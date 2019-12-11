/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Maize2000;

import analysis.maizeGeneticLoad.CrossMapUtils;
import utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author Aoyue
 */
public class Maize2000Entrance {
    public Maize2000Entrance(){
        this.otherClass();
        //this.mergeCountVariantFile();
       // this.crossMapConvert();
       
        
       
        
    }
    public void crossMapConvert(){
        /*如何利用public plant genetics中的 utils CrossMapUtils转换V3的坐标？？
        先将输入文件bed建立一个list,采用paralleStream方法一个一个多线程转换。
        */
        String myPythonPath = "/Users/Aoyue/miniconda3/bin/python";
        String myCrossMapPath = "/Users/Aoyue/miniconda3/bin/CrossMap.py";
        String myMaizeChainPath = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/AGPv3_to_AGPv4.chain.gz";
        String infileS = "/Users/Aoyue/Documents/1.bed";
        String outfileS = "/Users/Aoyue/Documents/1_V4.bed";
        CrossMapUtils cm = new CrossMapUtils(myPythonPath, myCrossMapPath, myMaizeChainPath,infileS, outfileS);
        cm.convert();
        
  
    }
    
    public void mergeCountVariantFile(){
        String infileDirS = "/Users/Aoyue/Documents/SIFTCalculate/ffffffff";
        String outfileS = "/Users/Aoyue/Documents/SIFTCalculate/allSIFTpredictions.txt";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        try{                       
            StringBuilder sb = new StringBuilder();
            sb.append("CHROM").append("\t").append("NONSynonymous").append("\t").append("SIFT_deleterious_mutation")
                        .append("\t").append("Non-synonymous_tolerent_mutation")
                        .append("\t").append("Neutral_mutation").append("\t").append("Startlost_mutation").append("\n");
            for(int i = 0; i < fs.length; i++){
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String subFileS = fs[i].getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(subFileS);
                String temp = null;
                while((temp = br.readLine()) != null){
                    if(temp.startsWith("CHROM")) continue;
                    sb.append(temp).append("\n");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();  
            }
        }
        catch(Exception e){
            System.out.println();
            e.printStackTrace();
            System.exit(1);              
        }
    }
    
   
        
    public void otherClass(){
        //new HapMap3Processorcp();
        //new HapMapTaxaProcessorcp();
        //new GERPcp();
        //new SIFTcp();
        new FastqQualitycp();
        
    }
    public static void main (String[] args) throws IOException{
        System.out.println("hello");
        new Maize2000Entrance();
     
        
       
    }
    
}

