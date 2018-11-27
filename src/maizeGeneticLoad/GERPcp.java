

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package maizeGeneticLoad;

import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import format.table.RowTable;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.CrossMapUtils;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author feilu
 */
class GERPcp {
    
    public GERPcp () {
        this.convertToV4GerpFile();
    }
    
    public void convertToV4GerpFile () {
        //this.mkAGPV3BED();
        //this.crossMapConvert();
        this.makeV3V4Map();
        this.mkV3V4mapBYBIXIAOYUE();
//        this.mkV4Gerp();
    }
    
     public void mkV4Gerp () {
        String inGerpDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/gerp/agpV3";
        String v3v4MapDirS = "/Users/feilu/Documents/database/maize/crossMap/V3V4Map";
        String outGerpDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/gerp/agpV4";
        File[] fs = new File(v3v4MapDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "map");
        List<File> fList = Arrays.asList(fs);
        HashIntIntMap[] v3v4Maps = new HashIntIntMap[fList.size()];
        fList.parallelStream().forEach(f -> {
            int index = Integer.valueOf(f.getName().replaceFirst(".map", "").replaceFirst("chr", ""))-1;
            v3v4Maps[index] = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
            int cnt = 0;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                    l = PStringUtils.fastSplit(temp);
                    if (l.get(2).startsWith("N")) {
                        v3v4Maps[index].put(Integer.parseInt(l.get(1)), -1);
                    }
                    else {
                        v3v4Maps[index].put(Integer.parseInt(l.get(1)), Integer.parseInt(l.get(2)));
                    }
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        fs = new File(inGerpDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "full");
        fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int index = Integer.valueOf(f.getName().replaceFirst(".msa.in.rates.full", "").replaceFirst("roast.chrom.", ""))-1;
            HashIntIntMap cMap = v3v4Maps[index];
            String outfileS = new File (outGerpDirS, "chr"+PStringUtils.getNDigitNumber(3, index+1)+"_Gerp.txt").getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Chr\tPos\tTreeLength\tValue");
                bw.newLine();
                int cnt = 0;
                String temp = null;
                List<String> l = null;
                int v4Pos = 0;
                StringBuilder sb = null;
                int chr = index + 1;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                    v4Pos = cMap.get(cnt);
                    if (v4Pos == -1) continue;
                    l = PStringUtils.fastSplit(temp);
                    sb = new StringBuilder();
                    sb.append(chr).append("\t").append(v4Pos).append("\t").append(l.get(0)).append("\t").append(l.get(1));
                    bw.write(sb.toString());
                    bw.newLine();
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
     
    public void mkV3V4mapBYBIXIAOYUE(){
        String V3infileDirs = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/temp1";
        String V4infileDirs = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/temp2";
        String outfileDirs = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/V3V4map";
        File[] f1 = new File(V4infileDirs).listFiles();
        f1 = IOUtils.listFilesEndsWith(f1, ".bed");
        List<File> fList = Arrays.asList(f1);
        fList.parallelStream().forEach(f -> {
            String V3infileS = new File(V3infileDirs, f.getName()).getAbsolutePath();
            String unmap = new File(V4infileDirs,f.getName().replaceFirst(".bed", ".bed.unmap")).getName();
            String outfileS = new File(outfileDirs, f.getName().replaceFirst(".bed.unmap", "_xiaoyue.map")).getAbsolutePath();
            //File[] unmap = IOUtils.listFilesEndsWith(new File(V4infileDirs).listFiles(),"unmap");
            //为了找到对应的unmap文件
            //List<File> unmapList = Arrays.asList(unmap);
            //int Index = Integer.parseInt(f.getName().replaceFirst("chr", "").replaceFirst(".map", ""));     
            List unmapList = null;
            String temp = null;
            try{
                BufferedReader br = IOUtils.getTextReader(unmap);
                while((temp = br.readLine()) != null){
                    temp = PStringUtils.fastSplit(temp).get(1);
                    unmapList.add(temp);
                }
                BufferedReader br1 = IOUtils.getTextReader(V3infileS);
                BufferedReader br2 = IOUtils.getTextReader(f.getName());
                while((temp = br1.readLine()) != null){
                    
                }
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
        });
    }
    
    public void makeV3V4Map () {
        String inDirS = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/temp1";
        String inCMDirS = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/temp2";
        String outDirS = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/V3V4map";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "bed");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {
            String outfileS = new File(outDirS, f.getName().replaceFirst(".bed", "_aoyue2Method.map")).getAbsolutePath();
            String unmapFileS = new File(inCMDirS, f.getName()+".unmap").getAbsolutePath();
            String mapFileS = new File(inCMDirS, f.getName()).getAbsolutePath();
            TIntArrayList unPosList = new TIntArrayList();
            
            try {
                /*先读入没有比对上的文件chr001.bed.unmap，将第2列的位置信息加入unPosList中，并转化为数组，排序；
                将chr001.bed V3版本的读入br，将chr001.bed V4版本的读入br2，将map文件写出 V3V4Map路径中，格式为 chr Pos_V3 Pos_V4,
                将V3版本的bed读入，第一列为染色体号，写入。第二列为V3的坐标系，赋值给query，并输出，
                将query在未比对成功的整型数组int[] unPos中查找，将二分搜索的值设置为index，
                如果找到，说明V3没有比对上，inedx>0, 在Pos_V4 第3列中输出NA，表示没有比对上；
                如果没有找到，说明V3比对上了，index<0, 进入循环 
                while (cnt < query) {
                            cnt++;
                            temp1 = br2.readLine();
                            l = PStringUtils.fastSplit(temp1);
                        }
                        sb.append(l.get(1));
                这个循环的解释：
                实践证明，不需要加循环。
                */
                BufferedReader br = IOUtils.getTextReader(unmapFileS); 
                String temp = null;
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    unPosList.add(Integer.parseInt(l.get(1)));
                }
                br.close();
                int[] unPos = unPosList.toArray(); /*将unPosList转化为整形数组，并按升序排序*/
                Arrays.sort(unPos);
                br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedReader br2 = IOUtils.getTextReader(mapFileS);
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Chr_V3\tPos_V3\tChr_V4\tPos_V4");
                bw.newLine();
                String temp1 = null;
                StringBuilder sb = new StringBuilder();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);  /*Temp是V3格式的文件坐标*/
                    sb = new StringBuilder(l.get(0));
                    int query = Integer.parseInt(l.get(1));
                    int index = Arrays.binarySearch(unPos, query);
                    int pos = query+1;
                    sb.append("\t").append(pos).append("\t");                   
                    if (index < 0)  {
                        //while (cnt < pos) {
                           // cnt++;
                            temp1 = br2.readLine();
                            l = PStringUtils.fastSplit(temp1);
                       // }
                        sb.append(l.get(0)).append("\t").append(Integer.parseInt(l.get(1))+1);
                    
                    }
                    else {
                        sb.append("NA").append("\t").append("NA");
                       // cnt++;
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                br2.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                //System.out.println(temp);
                e.printStackTrace();
            }
        });
    }
    
    public void crossMapConvert () {
        /*如何利用public plant genetics中的 utils CrossMapUtils转换V3的坐标？？
        先将输入文件bed建立一个list,采用paralleStream方法一个一个多线程转换。
        */
        String myPythonPath = "/Users/Aoyue/miniconda3/bin/python";
        String myCrossMapPath = "/Users/Aoyue/miniconda3/bin/CrossMap.py";
        String myMaizeChainPath = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/AGPv3_to_AGPv4.chain.gz";
        String intDirS = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/temp1";
        String outDirS = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/temp2";
        new File (outDirS).mkdir();
        File[] fs = new File(intDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "bed");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = new File(outDirS, f.getName()).getAbsolutePath();
            CrossMapUtils cm = new CrossMapUtils(myPythonPath, myCrossMapPath, myMaizeChainPath,f.getAbsolutePath(), outfileS);
            cm.convert();
        });
    }
        
   public void mkAGPV3BED () {
        String infileS = "/data1/home/aoyue/position/ChrLenCentPosi_agpV3.txt";
        String outDirS = "/data1/home/aoyue/position/temp";
        /*研究思路：表格读入并进去循环 -- 将第一列的染色体信息写进list整型类型中 -- 将第二列的染色体长度信息（一共10条）写入一个整型数组中
        -- 通过list我们可以利用parallelStream方法
        由于*/
        RowTable<String> t = new RowTable(infileS);
        List<Integer> chrList = new ArrayList();
        int[] chrLength = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            String s = t.getCell(i, 0);
            chrList.add(Integer.parseInt(s));
            chrLength[i] = Integer.parseInt(t.getCell(i, 1));
        }
        chrList.parallelStream().forEach(chr -> {
            String outfileS = new File (outDirS, "chr"+PStringUtils.getNDigitNumber(3, chr)+".bed").getAbsolutePath();
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < chrLength[chr-1]; i++) {
                    sb = new StringBuilder();
                    sb.append(chr).append("\t").append(i).append("\t").append(i);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
}
