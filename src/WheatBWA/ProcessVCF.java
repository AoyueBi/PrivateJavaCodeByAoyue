/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatBWA;

import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class ProcessVCF {
    public ProcessVCF (){
        this.sampleVCFdataInOrder();
        //this.sampleVCFdataRandom();
        //this.sampleVCFdataRandom2(); //该方法弃用
        
    }
    
    /**************** 根据库文件过滤后的位点数记录，在VCF文件中随机抽取以上记录的 int行，按顺序随机抽取**************/
    private void sampleVCFdataRandom2 () {
        //思想：随机抽取1万个位点
        /**hmpInfoFileS file
         * Chr	Pos	Ref	Alt	Ancestral(Gerp)	Major	Minor	MajorAlleleFrequency	MinorAlleleFrequency	SiteDepth	HetCount
            10	228730	A	T	NA	A	T	0.986692	0.013307985	3591	1
            10	228743	T	<DEL>	NA	T	<DEL>	0.9964455	0.0035545023	4072	0
            10	228744	A	<DEL>	NA	A	<DEL>	0.9964413	0.0035587188	4054	0
         */
        String hmpFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/hmp321_agpv4_chr10.vcf.gz";
        String hmpInfoFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/annoDB_new/hmp321Info_chr010_AGPv4_AnnoDB.txt.gz";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/hmp321_agpv4_chr10_10000.vcf";
        int snpNum = 0;
        int size = 10000;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(hmpInfoFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            int cnt = 0; 
            while ((temp = br.readLine()) != null) {
                cnt++; //第一行是2
            }
            snpNum = cnt;
            /*建立一个50000个 indices的整型数组， 赋值为indice[0]= */
            int[] indices = new int[size]; 
            for (int i = 0; i < size; i++) {
                indices[i] = (int)(Math.random()*snpNum); //随机挑选50000个数目
            }
            Arrays.sort(indices);
            br = IOUtils.getTextGzipReader(hmpFileS);
            cnt = 0;
            while ((temp = br.readLine()) != null) {
                if((temp = br.readLine()).startsWith("#")) {
                    bw.write(temp);bw.newLine();
                }
                cnt++; // 对snp开始计数
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" lines");
                int idx = Arrays.binarySearch(indices, cnt-1);
                if (idx < 0) continue;
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();bw.close();br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**************** 随机抽取VCF文件行数的 0.1%，按顺序随机抽取**************/
    private void sampleVCFdataRandom () {
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/hmp321_agpv4_chr10.vcf.gz";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/chr10_10.vcf";
        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if(temp.startsWith("#")) {
                    bw.write(br.readLine());
                    bw.newLine();
                }
                double r = Math.random();
                if (r > 0.0001) continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                List<String> l = PStringUtils.fastSplit(temp);
                if (l.get(3).contains(",")) continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void sampleVCFdataInOrder2 () {
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/hmp321_agpv4_chr10.vcf.gz";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/chr10_1000.vcf";
        int n = 1000;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            for (int i = 0; i < n; i++) {
                bw.write(br.readLine());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**************** 按顺序抽取VCF文件的前 int=lines 行，按顺序抽取，并且包含元文件**************/
    public void sampleVCFdataInOrder(){
        String infileS = "/Users/Aoyue/project/maizeGeneticLoad/oriData/hmp321_agpv4_chr10.vcf.gz";
        //String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/singleChrom/chr1_5000.vcf.gz";
        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/oriData/hmp321_agpv4/hmp321_agpv4_chr10.vcf.gz";
        int lines = 1000;
        int cnt = 0;
        try{
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            while((temp = br.readLine()) != null){
                cnt++;
                bw.write(temp);
                bw.newLine();
                if(cnt > lines) break;
            }
            bw.flush();bw.close();br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
}
