/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import static cern.clhep.Units.s;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class VCFtools {
    
    public VCFtools(){
        //this.calSNPHeter();
        this.calSNPMaf();
        
    }
    public VCFtools(String infileDirS, String outfileDirS){
        //this.calIndiHeter(infileDirS,outfileDirS);
    }
    
    public void calSNPMaf(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/000_sampleVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/004_calMaf";
        File[] fs = new File(infileDirS).listFiles();
        //fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            //开始处理文件，计时
            long startTime = System.nanoTime();
            Calendar cal = Calendar.getInstance();
            int hour = cal.get(Calendar.HOUR_OF_DAY);
            int minute = cal.get(Calendar.MINUTE);
            int second = cal.get(Calendar.SECOND);
            System.out.println("******************************************************" );
            System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + "Now starting to " + f.getName() + " minor allele frequency;"); 
            
            try{
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS,f.getName().split(".vc")[0] + "_maf.txt").getAbsolutePath();
               
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")){
                    br = IOUtils.getTextReader(infileS);
                }
                else if (infileS.endsWith(".vcf.gz")){
                    br = IOUtils.getTextGzipReader(infileS);
                }
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                //bw.write("Chr\tPos\tHetPropotion\n");
                bw.write("Chr\tPos\tMinorAlleleGametes\tMaf\n");
                String temp;
                String te[] = null;
                while ((temp=br.readLine())!=null){
                    double genoNum = 0;
                    double refAlleleGametes = 0;
                    double altAlleleGametes = 0;
                    double refAF = 0;
                    double altAF = 0;
                    double maf = 0;
                    //在一个位点内进行计算
                    if(!temp.startsWith("#")){
                        te = temp.split("\t");
                        if(te[4].length() >1){
                            continue;
                        }
                        for (int i = 9; i < te.length;i++){
                            if(!te[i].startsWith(".")){
                                genoNum++; //have the genotype
                                if(te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
                                    refAlleleGametes++; 
                                    altAlleleGametes++; 
                                }
                                if(te[i].startsWith("0/0")) {
                                    refAlleleGametes++;
                                    refAlleleGametes++;
                                }
                                if(te[i].startsWith("1/1")) {
                                    altAlleleGametes++;
                                    altAlleleGametes++;
                                }
                            }
                        }
                        refAF = refAlleleGametes/(refAlleleGametes+altAlleleGametes);
                        altAF = altAlleleGametes/(refAlleleGametes+altAlleleGametes);;
                        if(refAF >= altAF){
                                maf = altAF;
                        }else {
                                maf = refAF;
                        }
                        //bw.write(te[0]+ "\t"+ te[1]+ "\t" + String.format("%.5f", hetRate) + "\n");
                        bw.write(te[0]+ "\t"+ te[1]+ "\t" + String.format("%.0f", altAlleleGametes)+ "\t"+ String.format("%.5f", maf) +"\n");
                    }
                }
                br.close();bw.flush();bw.close();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
            
            //文件处理完毕，计时
                hour = cal.get(Calendar.HOUR_OF_DAY);
                minute = cal.get(Calendar.MINUTE);
                second = cal.get(Calendar.SECOND);
                long endTime = System.nanoTime();
                float excTime = (float) (endTime - startTime) / 1000000000;
                //System.out.println("******************************************************" );
                //System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + f.getName() + " is finished!!!");
                System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
        });
    }
    
    /**
     * Calculate the heterozygote count and propotion by snp site
     */
    public void calSNPHeter(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/000_sampleVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/002_calSNPHeter";
        File[] fs = new File(infileDirS).listFiles();
        //fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            //开始处理文件，计时
            long startTime = System.nanoTime();
            Calendar cal = Calendar.getInstance();
            int hour = cal.get(Calendar.HOUR_OF_DAY);
            int minute = cal.get(Calendar.MINUTE);
            int second = cal.get(Calendar.SECOND);
            System.out.println("******************************************************" );
            System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + "Now starting to " + f.getName() + " heterozygote propotion;"); 
            
            try{
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS,f.getName().split(".vc")[0] + "_SNPheter.txt").getAbsolutePath();
               
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")){
                    br = IOUtils.getTextReader(infileS);
                }
                else if (infileS.endsWith(".vcf.gz")){
                    br = IOUtils.getTextGzipReader(infileS);
                }
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                //bw.write("Chr\tPos\tHetPropotion\n");
                bw.write("Chr\tPos\tHetNum\tHomNum\tHetPropotion\n");
                String temp;
                String te[] = null;
                while ((temp=br.readLine())!=null){
                    int genoNum = 0;
                    double homNum = 0;
                    double hetNum = 0;
                    double hetRate = 0;
                    //在一个位点内进行计算
                    if(!temp.startsWith("#")){
                        te = temp.split("\t");
                        if(te[4].length() >1){
                            continue;
                        }
                        for (int i = 9; i < te.length;i++){
                            if(!te[i].startsWith(".")){
                                genoNum++; //have the genotype
                                if(te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
                                    hetNum++; //the number of heterozygous
                                }
                                if(te[i].startsWith("0/0") || te[i].startsWith("1/1")) {
                                    homNum++; //the number of heterozygous
                                }
                            }
                        }
                        hetRate = hetNum/genoNum;
                        //bw.write(te[0]+ "\t"+ te[1]+ "\t" + String.format("%.5f", hetRate) + "\n");
                        bw.write(te[0]+ "\t"+ te[1]+ "\t" + String.format("%.0f", hetNum)+ "\t"+ String.format("%.0f", homNum) + "\t"+ String.format("%.5f", hetRate) + "\n");
                    }
                }
                br.close();bw.flush();bw.close();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
            
            //文件处理完毕，计时
                hour = cal.get(Calendar.HOUR_OF_DAY);
                minute = cal.get(Calendar.MINUTE);
                second = cal.get(Calendar.SECOND);
                long endTime = System.nanoTime();
                float excTime = (float) (endTime - startTime) / 1000000000;
                //System.out.println("******************************************************" );
                //System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + f.getName() + " is finished!!!");
                System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
        });
    }
    
    /**
     * Calculate the heterozygote count and propotion by individual taxa
     */
    public void calIndiHeter(String infileDirS,String outfileDirS){
        //String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/000_sampleVCF";
        //String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/002_calSNPHeter";
        File[] fs = new File(infileDirS).listFiles();
        //fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            //Start to cal the time beginning
            long startTime = System.nanoTime();
            Calendar cal = Calendar.getInstance();
            int hour = cal.get(Calendar.HOUR_OF_DAY);
            int minute = cal.get(Calendar.MINUTE);
            int second = cal.get(Calendar.SECOND);
            System.out.println("******************************************************" );
            System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + "Now starting to " + f.getName() + " heterozygote propotion;"); 
            //Start to cal heterozygous
            try{
                String infileS = f.getAbsolutePath();
                String siteoutfileS = new File(outfileDirS,f.getName().split(".vc")[0] + "_SNPheter.txt").getAbsolutePath();
                String indioutfileS = new File(outfileDirS,f.getName().split(".vc")[0] + "_IndiHeter.txt").getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".vcf"))br = IOUtils.getTextReader(infileS);
                else if (infileS.endsWith(".vcf.gz"))br = IOUtils.getTextGzipReader(infileS);
                BufferedWriter bw = IOUtils.getTextWriter(siteoutfileS);
                BufferedWriter indibw = IOUtils.getTextWriter(indioutfileS);
                bw.write("Chr\tPos\tHetNum\tHomNum\tHetPropotion\n");
                
                String temp;
                String te[] = null;
                /*******************定义taxa数组，并添加元素***********************/
                while((temp=br.readLine()).startsWith("##")) {}
                List<String> l = PStringUtils.fastSplit(temp);
                String[] taxa = new String[l.size()-9];
                for(int i =0; i<taxa.length;i++){
                    taxa[i]= l.get(i+9);
                }
                /*******************定义Genotype TIntArrayList***********************/
                //TIntArrayList[] genoList = new TIntArrayList[taxa.length];
                //for (int i = 0; i < taxa.length; i++) genoList[i] = new TIntArrayList();
                List[] genoList = new ArrayList[taxa.length];
                for (int i = 0; i < taxa.length; i++) genoList[i] = new ArrayList();
                
 
                while ((temp=br.readLine())!=null){
                    int genoNum = 0;
                    double homNum = 0;
                    double hetNum = 0;
                    double hetRate = 0;
                    //在一个位点内进行计算
                    
                    if(!temp.startsWith("#")){
                        //l = PStringUtils.fastSplit(temp, "\t");
                        te = temp.split("\t");
                        if(te[4].length() >1){
                            continue;
                        }
                        for (int i = 9; i < te.length;i++){
                            if(te[i].startsWith(".")){genoList[i-9].add(0);}
                            if(!te[i].startsWith(".")){
                                genoNum++; //have the genotype
                                if(te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
                                    hetNum++; //the number of heterozygous
                                    genoList[i-9].add(2); // 2 stand for heter
                                }
                                if(te[i].startsWith("0/0") || te[i].startsWith("1/1")) {
                                    homNum++; //the number of heterozygous
                                    genoList[i-9].add(1); // 1 stand for homo
                                }
                            }
                        }
                        hetRate = hetNum/genoNum;
                        bw.write(te[0]+ "\t"+ te[1]+ "\t" + String.format("%.0f", hetNum)+ "\t"+ String.format("%.0f", homNum) + "\t"+ String.format("%.5f", hetRate) + "\n");
                    }
                }
               
                
                /******************* 开始计算个体的杂合度***********************/
                indibw.write("INDV\tTotalSitesWithGeno\tHetSites\tHetProportion\tMissingSites\tMissProportion\n");
                for (int i = 0; i < taxa.length; i++) {
                    double hetSites =  Collections.frequency(genoList[i], 2);
                    double totalSitesWithGeno =  Collections.frequency(genoList[i], 1) + Collections.frequency(genoList[i], 2); //含有基因型的位点数
                    double hetProportion = hetSites/totalSitesWithGeno;
                    double missSites = Collections.frequency(genoList[i], 0);
                    double missProportion = missSites/(Collections.frequency(genoList[i], 0)+ Collections.frequency(genoList[i], 1) + Collections.frequency(genoList[i], 2));
                    indibw.write(taxa[i]+"\t" + String.format("%.0f", totalSitesWithGeno) + "\t"+ String.format("%.0f", hetSites) + "\t"  + String.format("%.5f", hetProportion) + "\t" +
                            String.format("%.0f", missSites) + "\t" + String.format("%.5f", missProportion) + "\n");
                }
                br.close();bw.flush();bw.close();
                indibw.flush();
                indibw.close();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
            
            //文件处理完毕，计时
                hour = cal.get(Calendar.HOUR_OF_DAY);
                minute = cal.get(Calendar.MINUTE);
                second = cal.get(Calendar.SECOND);
                long endTime = System.nanoTime();
                float excTime = (float) (endTime - startTime) / 1000000000;
                //System.out.println("******************************************************" );
                //System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + f.getName() + " is finished!!!");
                System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
                
                /**
                 * ******************************************************
                    Here is the main class of Basic Genetics Analysis !
                    ******************************************************
                    22:54:0	Now starting to chr001_5000.vcf heterozygote propotion;
                    Execution time: 0.39s
                    ******************************************************
                    22:54:0	Now starting to chr006_5000.vcf heterozygote propotion;
                    Execution time: 0.29s
                    ******************************************************
                    22:54:1	Now starting to chr002_5000.vcf heterozygote propotion;
                    Execution time: 0.19s
                    ******************************************************
                    22:54:1	Now starting to chr003_5000.vcf heterozygote propotion;
                    Execution time: 0.13s
                    ******************************************************
                    22:54:1	Now starting to chr005_5000.vcf heterozygote propotion;
                    Execution time: 0.10s
                    ******************************************************
                    22:54:1	Now starting to chr004_5000.vcf heterozygote propotion;
                    Execution time: 0.10s
                 */
        });
    }
}
