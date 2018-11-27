

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package maizeGeneticLoad;

import format.genomeAnnotation.GeneFeature;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author feilu
 */
class SIFTcp {
    
    public SIFTcp () {
        //this.testData();
        this.addSIFT();
    }
    
    public void addSIFT () {
        String infileS = "/Users/Aoyue/Documents/Data/referenceGenome/GeneAnnotation/Zea_mays.AGPv4.38.pgf";
        String siftDirS = "/Users/Aoyue/Documents/SIFTCalculate/01SIFTresult/allresult";
        String annotationDBDirS = "/Volumes/LuLab3T_30/maizeLoad/annoDB";
        String outputDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/annoDB_newnew";
        new File (outputDirS).mkdir();
        GeneFeature gf = new GeneFeature(infileS);
        int geneNum = gf.getGeneNumber();
        String[] trans = new String[geneNum];
        for (int i = 0; i < geneNum; i++) {
            int index = gf.getLongestTranscriptIndex(i);
            trans[i] = gf.getTranscriptName(i, index); /*得到第一个基因，最长转录本的index的转录本名字*/
        }
        Arrays.sort(trans);
        File[] fs = new File(siftDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "txt");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().split("_")[2].replaceFirst("chr", ""))-1;
//            String infileSS = "/Users/feilu/Documents/analysisL/production/maizeLoad/sift/source/hmp321_agpv4_chr10_SIFTannotations.xls";
//            f = new File(infileSS);
//            chrIndex = 9;
            //hmp321Info_chr001_AGPv4_AnnoDB.txt
            String dbFileS = new File (annotationDBDirS, "hmp321Info_chr"+PStringUtils.getNDigitNumber(3, chrIndex+1)+"_AGPv4_AnnoDB.txt").getAbsolutePath();
            String outfileS = new File(outputDirS,  "hmp321Info_chr"+PStringUtils.getNDigitNumber(3, chrIndex+1)+"_AGPv4_AnnoDB.txt").getAbsolutePath();
            String header = null;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                List<SIFTRecord> sList = new ArrayList<>(); /*list也可以是自己创造的类型，不是非要特指Integer String 等这些类*/
                List<String> l = null;
                /*将文件读入，并把每一条记录放入 sList中*/
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String alt = l.get(3);
                    String transcript = l.get(4).replaceFirst("transcript.", "");
                    String type = l.get(8);
                    String value = l.get(12);
                    SIFTRecord s = new SIFTRecord( pos, alt,  transcript,  type, value);
                    sList.add(s);
                }
                br.close();
                Collections.sort(sList);/*要想进行搜索，被搜索的sList必须是按自然顺序排序的*/
                /*将数据库dbFileS中的文件读入，建立query*/
                br = IOUtils.getTextReader(dbFileS);
                header = br.readLine()+"\tVariant_type\tSIFT_score\tTranscript";
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write(header);
                bw.newLine();
                int index = -1;
                int startIndex = -1;
                int endIndex = -1;
                StringBuilder sb = null;
                while ((temp = br.readLine()) != null) {
                    sb = new StringBuilder();
                    sb.append(temp).append("\t");
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String alt = l.get(3).split(",")[0];
                    SIFTRecord query = new SIFTRecord (pos, alt, "", "", "");
                    index = Collections.binarySearch(sList, query);
                    if (index < 0) { /*如果小于0，说明位置和alt碱基不匹配，因为sift结果中有很多位置pos不存在，就算位置存在，alt碱基也可能和header库中的不一样，这些情况的site都没有sift值*/
                        sb.append("NA").append("\t").append("NA").append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    startIndex = index;
                    endIndex = index;
                    while ((startIndex-1)>-1 && sList.get(startIndex-1).isSimilar(pos, alt)) {
                        startIndex--;
                    }
                    while ((endIndex+1) < sList.size() && sList.get(endIndex+1).isSimilar(pos, alt)) {
                        endIndex++;
                    }
                    boolean status = false;
                    for (int i = startIndex; i < endIndex+1; i++) {
                        if (Arrays.binarySearch(trans, sList.get(i).transcript) >= 0 && sList.get(i).type.equals("NONSYNONYMOUS") || sList.get(i).type.equals("SYNONYMOUS")) {
                            if (Arrays.binarySearch(trans, sList.get(i).transcript) < 0) {
                                sb = new StringBuilder();
                                sb.append(temp).append("\t");sb.append("NA\tNA\tNA");
                            }
                            else {
                                sb = new StringBuilder();
                                sb.append(temp).append("\t");
                                sb.append(sList.get(i).type).append("\t").append(sList.get(i).value).append("\t").append(sList.get(i).transcript);
                            }
                            
                            bw.write(sb.toString());
                            bw.newLine();
                            status = true;
                            break;
                        }
                    }
                    if (status ==false) {
                        sb = new StringBuilder();
                        sb.append(temp).append("\t");
                        sb.append("NA\tNA\tNA");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    class SIFTRecord implements Comparable<SIFTRecord>{
        public int pos;
        public String alt;
        public String transcript;
        public String type;
        public String value;
        
        public SIFTRecord (int pos, String alt, String transcript, String type, String value) {            
            this.pos = pos;
            this.alt = alt;
            this.transcript = transcript;
            this.type = type;
            this.value = value;
        }
        
        public boolean isSimilar (int pos, String alt) {
            if (pos == this.pos && alt.equals(this.alt)) return true;
            return false;
        }
        
        @Override
        public int compareTo(SIFTRecord o) {
            if (this.pos < o.pos) {
                return -1;
            }
            else if (this.pos == o.pos) {
                return this.alt.compareTo(o.alt);
            }
            else {
                return 1;
            }
            
        }
    }
    
    public boolean ifMatch (List<String> l, int pos, String alt, String[] trans) {
        int posS = Integer.parseInt(l.get(1));
        if (posS != pos) return false;
        if (!alt.equals(l.get(3))) return false;
        String type = l.get(8);
        if (!(type.startsWith("NONSYNONYMOUS") || type.startsWith("SYNONYMOUS"))) return false;
        String query = l.get(4).replaceFirst("transcript.", "");
        if (Arrays.binarySearch(trans, query) < 0) return false;
        return true;
    }
    
    public void testData () {
        String infileS = "/Users/feilu/Documents/analysisL/production/maizeLoad/sift/source/hmp321_agpv4_chr10_SIFTannotations.xls";
        Set<String> typeSet = new HashSet<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                typeSet.add(tem[8]);
            }
            br.close();
            String[] types = typeSet.toArray(new String[typeSet.size()]);
            Arrays.sort(types);
            Set<String>[] subsets = new HashSet[types.length];
            for (int i = 0; i < subsets.length; i++) subsets[i] = new HashSet<>();
            br = IOUtils.getTextReader(infileS);
            temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                String query = tem[8];
                int index = Arrays.binarySearch(types, query);
                subsets[index].add(tem[16]);
            }
            br.close();
            for (int i = 0; i < types.length; i++) {
                System.out.println(types[i]);
                System.out.println(subsets[i]+"\n");
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
}
