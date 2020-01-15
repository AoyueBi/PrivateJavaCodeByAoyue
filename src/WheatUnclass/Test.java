/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatUnclass;

import format.table.RowTable;
import gnu.trove.list.array.TIntArrayList;
import utils.IOUtils;
import utils.PStringUtils;
import utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class Test {
    int a = 4;
    
    public Test (int b) {
        this.a = b;
        System.out.println(this.a);
        System.out.println("#####");
        System.out.println(a);
    }


    /**
     *
     */

    public void mkBurdenFiles () {
        //本块内容进行文件夹的建立
        String vcfDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/002_exonSNPVCF";
        String taxaFileS = "/Users/feilu/Documents/analysisH/vmap2/006_populationStructure/vmap2_taxa.txt";
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousISite";
        File[] fs = new File(dirS).listFiles();
        List<File> dirList = new ArrayList<>();
        for (int i = 0; i < fs.length; i++) {
            if (!fs[i].isDirectory()) continue;
            dirList.add(fs[i]);
        }
        Collections.sort(dirList);


        //本块内容进行分组的建立
        RowTable<String> t = new RowTable<>(taxaFileS);
        int groupIndex = t.getColumnIndex("TreeValidatedGroupbyPloidy"); //即倍性分组那一列所在的index
        List<String> gList = t.getColumn(groupIndex);
        HashSet<String> groupSet = new HashSet<>(gList);
        groupSet.remove("ExclusionHexaploid");
        groupSet.remove("ExclusionTetraploid");
        String[] groups = groupSet.toArray(new String[groupSet.size()]);
        Arrays.sort(groups);


        //this section 建立按照倍性分组的taxa列表，每个分组有个 list
        List<String>[] groupTaxaLists = new ArrayList[groups.length];
        String[] groupPloidy = new String[groups.length];
        for (int i = 0; i < groupTaxaLists.length; i++) {
            groupTaxaLists[i] = new ArrayList<>();
        }
        int index = -1;
        String query = null;
        for (int i = 0; i < t.getRowNumber(); i++) { //进入种质信息表中进行循环
            query = t.getCell(i, groupIndex); //得到分组所在列的分组信息
            index = Arrays.binarySearch(groups, query);
            if (index < 0) continue; //小于0说明是 ExclusionHexaploid or ExclusionTetraploid 不进行分析
            groupTaxaLists[index].add(t.getCell(i, 4)); //得到每个分组的taxa list
            groupPloidy[index] = t.getCell(i, 5); //得到倍性 AABBDD AABB DD
        }
        for (int i = 0; i < groupTaxaLists.length; i++) {
            Collections.sort(groupTaxaLists[i]);
        }
        dirList.parallelStream().forEach(dir -> {
            this.outputBurden(dir, vcfDirS, groups, groupTaxaLists, groupPloidy, dirS);
        });

    }


    /**
     * 解析代码
     * @param delDir
     * @param vcfDirS
     * @param groups
     * @param groupTaxaLists
     * @param groupPloidy
     * @param outputDirS
     */
    private void outputBurden (File delDir, String vcfDirS, String[] groups, List<String>[] groupTaxaLists, String[] groupPloidy, String outputDirS) {
        String outFileS = new File (outputDirS, delDir.getName()+".txt").getAbsolutePath();
        String[] genomeTypes = {"A", "B", "D"};

        try {
            BufferedWriter bw = IOUtils.getTextWriter(outFileS);
            String header = "Subgenome\tPloidy\tTaxa\tDerivedAlleleCount\tScoredGenotypeCount\tDerivedAllelePerSite";
            bw.write(header);
            bw.newLine();
            int[] groupIndex = new int[2];

            //对A B D 进行外层循环
            //对 AABBDD AABB DD进行第二小层循环
            //A --- 如果AABBDD 包含 A ，那么 groupIndex[0]=0; 如果 AABB包含A，那么groupIndex[1]=1；
            //B --- 如果AABBDD 包含 B ，那么 groupIndex[0]=0; 如果 AABB包含B，那么groupIndex[1]=1；
            //D --- 如果AABBDD 包含 D ，那么 groupIndex[0]=0; 如果 DD包含D，那么groupIndex[1]=1；

            for (int i = 0; i < genomeTypes.length; i++) { //A B D各自循环一遍
                int cnt = 0;
                for (int j = 0; j < groupPloidy.length; j++) {  //groupPloidy数组，包含AABBDD AABB DD
                    if (groupPloidy[j].contains(genomeTypes[i])) {
                        groupIndex[cnt] = j;
                        cnt++;
                    }
                }


                int[][] alleleCount = new int[groupIndex.length][];
                int[][] genotypeCount = new int[groupIndex.length][];
                int[][] taxaIndex = new int[groupIndex.length][];

                //对groupIndex[cnt]进行分析，如果groupIndex[0]，说明是AABBDD；如果groupIndex[1]，说明是AABB
                //AABBDD的 alleleCount[][]是第
                //外层循环A亚基因组，j=0,groupIndex[0]是0，AABBDD组，groupTaxaLists[0]是AABBDD组的taxa集合 alleleCount[0]是AABBDD的所有taxa的大小，即每个taxa有多少个derived allele count
                //该段循环是对每个taxa进行属性鉴定，查看其有多少个 有害突变，多少个基因型，多少个 taxaIndex
                for (int j = 0; j < groupIndex.length; j++) {
                    alleleCount[j] = new int[groupTaxaLists[groupIndex[j]].size()];
                    genotypeCount[j] = new int[groupTaxaLists[groupIndex[j]].size()];
                    taxaIndex[j] = new int[groupTaxaLists[groupIndex[j]].size()];
                }
                int[] chrIDs = RefV1Utils.getChrIDsOfSubgenome(genomeTypes[i]);
                for (int j = 0; j < chrIDs.length; j++) {
                    String vcfFileS = "chr"+PStringUtils.getNDigitNumber(3, chrIDs[j])+"_exon_vmap2.1.vcf.gz";
                    vcfFileS = new File (vcfDirS, vcfFileS).getAbsolutePath();
                    String delFileS = "chr"+PStringUtils.getNDigitNumber(3, chrIDs[j])+"_SNP_anno.txt.gz";
                    delFileS = new File (delDir, delFileS).getAbsolutePath();
                    String temp = null;
                    List<String> l = new ArrayList<>();
                    TIntArrayList posList = new TIntArrayList();
                    List<String> ancestraList = new ArrayList<>();
                    BufferedReader br = IOUtils.getTextGzipReader(delFileS); //只读取有害突变的位点和祖先状态
                    br.readLine();
                    String ances = null;
                    while((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp);
                        ances = l.get(15);
                        if (ances.startsWith("N")) continue;
                        posList.add(Integer.parseInt(l.get(2)));
                        ancestraList.add(ances);
                    }
                    br.close();
                    int[] delPos = posList.toArray();
                    String[] ancestral = ancestraList.toArray(new String[ancestraList.size()]);
                    br = IOUtils.getTextGzipReader(vcfFileS);
                    while ((temp = br.readLine()).startsWith("##")) {}
                    l = PStringUtils.fastSplit(temp);
                    int index = -1;
                    for (int k = 0; k < groupIndex.length; k++) {
                        for (int m = 0; m < l.size(); m++) {
                            index = Collections.binarySearch(groupTaxaLists[groupIndex[k]], l.get(m));
                            if (index < 0) continue;
                            taxaIndex[k][index] = m;
                        }
                    }
                    int derivedCode = -1;
                    String currentVCF = null;
                    String[] tem = null;
                    while ((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp.substring(0,60));
                        index = Arrays.binarySearch(delPos, Integer.parseInt(l.get(1)));
                        if (index < 0) continue;
                        l = PStringUtils.fastSplit(temp);
                        if (l.get(3).equals(ancestral[index])) derivedCode = 1;
                        else derivedCode = 0;
                        for (int k = 0; k < groupIndex.length; k++) {
                            for (int m = 0; m < groupTaxaLists[groupIndex[k]].size(); m++) {
                                currentVCF = l.get(taxaIndex[k][m]);
                                if (currentVCF.startsWith(".")) {
                                    continue;
                                }
                                genotypeCount[k][m]++;
                                tem = currentVCF.split(":")[0].split("/");
                                for (int n = 0; n < tem.length; n++) {
                                    if (Integer.parseInt(tem[n]) == derivedCode) {
                                        alleleCount[k][m]++;
                                    }
                                }
                            }
                        }
                    }
                    br.close();
                    System.out.println(vcfFileS);
                }
                StringBuilder sb = new StringBuilder();
                for (int j = 0; j < groupIndex.length; j++) {
                    for (int k = 0; k < alleleCount[j].length; k++) {
                        sb.setLength(0);
                        sb.append(genomeTypes[i]).append("\t").append(groupPloidy[groupIndex[j]]).append("\t");
                        sb.append(groupTaxaLists[groupIndex[j]].get(k)).append("\t");
                        sb.append(alleleCount[j][k]).append("\t").append(genotypeCount[j][k]).append("\t");
                        sb.append((float)((double)alleleCount[j][k]/genotypeCount[j][k]));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            bw.close();

        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void main (String[] args) {
        new Test(7);
    }
    
}