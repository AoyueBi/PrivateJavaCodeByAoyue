/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGBS;

import xujun.analysis.rnaseq.*;
import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import utils.IOUtils;
import xuebo.analysis.annotation.FStringUtils;

/**
 *
 * @author xujun
 */
public class GBScp {
    public GBScp(){
        this.barcode();
    }
    private void barcode () {
        String sampleInfor="/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/02IlluminaSeq/20180601_GBSLibrarybarcode.txt";
        RowTable rt=new RowTable(sampleInfor);
        HashMap hm1=new HashMap();
        HashMap hm2=new HashMap();
        List sample=new ArrayList();
        for(int i=0;i<rt.getRowNumber();i++){
            hm1.put(rt.getCell(i, 2), rt.getCell(i, 3));
            hm2.put(rt.getCell(i, 4), rt.getCell(i, 2));
            sample.add(rt.getCell(i, 2)); /*Plate编号*/
        }
        String infileDirS = "/Users/Aoyue/Downloads/2.cleandata/20180601-51-NEB12_TKD180600155";
        String outfileDirS = "/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/02IlluminaSeq/result";
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        int[] cnt=new int[sample.size()];
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_TKD180600155_1.clean.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_TKD180600155_1.clean.fq.gz").getAbsolutePath();
            String outfile = new File (outfileDirS, "count.txt").getAbsolutePath();
            int count=0;int count1=0;
            String temp1=null;String seq1=null;
            String temp2=null;String seq2=null;           
            try {
                BufferedReader br1 = utils.IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = utils.IOUtils.getTextGzipReader(infile2);   
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                while ((temp1 = br2.readLine()) != null) {
                    count1++;
                    seq1=br2.readLine(); /*读取FASTQ文件的第二行，存为序列 seq1*/
                    if(hm2.get(seq1.substring(0,6))==null){ /*hm2的key为MSP1的序列，value为Plate编号。通过序列的前6个碱基，找到Plate编号。  如果找到编号，则条件为假，如果找不到编号，则条件为真*/
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                        count++;
                    }else{
                        temp2=br1.readLine(); /*文件1的第一行信息，不是序列*/
                        seq2=br1.readLine(); /*文件1的第二行信息，是序列*/
                        if(hm1.get(hm2.get(seq1.substring(0,6)))!=null){ /*hm2.get(seq1.substring(0,6)不为空，能找到对应的Plate编号， 根据编号再找到相应的BamHI_barcode，如果
                            BamHI_barcode不为空*/
                            int index=sample.indexOf(hm2.get(seq1.substring(0,6))); /*找到改编号的索引*/
                            cnt[index]++; /**/
                        }else{
                            count++;
                        }
                        br1.readLine();br1.readLine();
                    }
                    br2.readLine();br2.readLine();
                }
                for(int i=0;i<cnt.length;i++){
                    bw.write(sample.get(i)+"\t"+cnt[i]);
                    bw.newLine();
                }
                System.out.println(count);
                System.out.println(count1);
                br1.close();br2.close();
                bw.flush();bw.close();
            }
            
            catch (Exception e) {
                e.printStackTrace();
            }

        });
    }
}
