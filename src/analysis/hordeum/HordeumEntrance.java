/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.hordeum;

import aoyue.analysis.MaizeGeneticLoad.Rangescp;
import format.range.Range;
import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class HordeumEntrance {
    
    public HordeumEntrance(){
        //this.divideIfGene();
        this.divideIfExon();
        
    }
    
    public void divideIfExon(){
        String infileS = "/Users/Aoyue/Downloads/barly/000_WithinGeneRNA.txt";
        String gff3FileS = "/Users/Aoyue/Downloads/barly/Hordeum_vulgare.IBSC_v2.42.gff3.gz";
        String withinExonFileS = "/Users/Aoyue/Downloads/barly/001_WithinExonRNA.txt";
        String withinIntronFileS = "/Users/Aoyue/Downloads/barly/001_withinIntronRNA.txt";
        RowTable<String> t = new RowTable<>(infileS);
        int cnt = 0;
        int cntwithinExon = 0;
        int cntwithinIntron =0;
        try{
            BufferedReader br = IOUtils.getTextGzipReader(gff3FileS);
            BufferedWriter bw = IOUtils.getTextWriter(withinExonFileS);
            BufferedWriter bw1 = IOUtils.getTextWriter(withinIntronFileS);
            bw.write("Chr\tStart\tEnd\tTotalCount");bw.newLine();
            bw1.write("Chr\tStart\tEnd\tTotalCount");bw1.newLine();
            String temp = br.readLine();
            List<String> l = null;
            ArrayList<Range> rangesList = new ArrayList();
            while((temp = br.readLine()) != null){
                if(temp.startsWith("#")) continue;
                if(temp.startsWith("Mt")) continue;
                if(temp.startsWith("Pt")) continue;
                if(temp.startsWith("chrUn")) continue;
                l = PStringUtils.fastSplit(temp);
                int chrdb = Integer.parseInt(l.get(0).split("r")[1].split("H")[0]);
                String biotype = l.get(2);
                int startdb = Integer.parseInt(l.get(3));
                int enddb = Integer.parseInt(l.get(4));
                if(biotype.equals("mRNA")) continue;
                if(biotype.equals("gene")) continue;
                if(biotype.equals("three_prime_UTR")) continue;
                if(biotype.equals("five_prime_UTR")) continue;
                if(biotype.equals("CDS")) continue;
                Range e = new Range(chrdb,startdb,enddb);
                rangesList.add(e);
            }
         
            System.out.println("The number of exon ranges is " +  rangesList.size()  );
            
            Rangescp db = new Rangescp(rangesList,"exon");
            db.sortByStartPosition();
            
        for(int i = 0; i<t.getRowNumber(); i++){
        String chrs = t.getCellAsString(i, 0);
        //if(chrs.startsWith("bgh")) continue;
        //if(chrs.startsWith("chrUn")) continue;
        cnt++;
        //chrs = chrs.split("r")[1].split("H")[0];
        int chr = Integer.parseInt(chrs);
        String starts = t.getCellAsString(i, 1); int start = Integer.parseInt(starts);
        String ends = t.getCellAsString(i, 2); int end = Integer.parseInt(ends);
        //String count = t.getCellAsString(i, 3).split("_")[1];
        String count = t.getCellAsString(i, 3);

        if((db.isInRanges(chr,start)) || (db.isInRanges(chr,end))){
            cntwithinExon++;
            bw.write(chrs + "\t" + starts + "\t" + ends + "\t" + count);
            bw.newLine();
        }
        else{
            cntwithinIntron++;
            bw1.write(chrs + "\t" + starts + "\t" + ends + "\t" + count);
            bw1.newLine();
        }
        }
        System.out.println("The number of reads within chr1H to chr7H gene region  is " +  cnt );
        System.out.println("The number of withinExon reads  is " +  cntwithinExon );
        System.out.println("The number of withinIntron reads is " +  cntwithinIntron );
        br.close();bw.flush();bw.close();bw1.flush();bw1.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void divideIfGene(){
        String infileS = "/Users/Aoyue/Downloads/barly/AllbarleyandBghsRNA.txt";
        String gff3FileS = "/Users/Aoyue/Downloads/barly/Hordeum_vulgare.IBSC_v2.42.gff3.gz";
        String withingeneFileS = "/Users/Aoyue/Downloads/barly/000_WithinGeneRNA.txt";
        String intergenicFileS = "/Users/Aoyue/Downloads/barly/000_intergenicRNA.txt";
        RowTable<String> t = new RowTable<>(infileS);
        int cnt = 0;
        int cntwithingene = 0;
        int cntintergenic =0;
        try{
            BufferedReader br = IOUtils.getTextGzipReader(gff3FileS);
            BufferedWriter bw = IOUtils.getTextWriter(withingeneFileS);
            BufferedWriter bw1 = IOUtils.getTextWriter(intergenicFileS);
            bw.write("Chr\tStart\tEnd\tTotalCount");bw.newLine();
            bw1.write("Chr\tStart\tEnd\tTotalCount");bw1.newLine();
            String temp = br.readLine();
            List<String> l = null;
            ArrayList<Range> rangesList = new ArrayList();
            while((temp = br.readLine()) != null){
                if(temp.startsWith("#")) continue;
                if(temp.startsWith("Mt")) continue;
                if(temp.startsWith("Pt")) continue;
                if(temp.startsWith("chrUn")) continue;
                //String tem = temp.substring(0, 91);
                l = PStringUtils.fastSplit(temp);
                int chrdb = Integer.parseInt(l.get(0).split("r")[1].split("H")[0]);
                String biotype = l.get(2);
                int startdb = Integer.parseInt(l.get(3));
                int enddb = Integer.parseInt(l.get(4));
                if(biotype.equals("mRNA")) continue;
                if(biotype.equals("exon")) continue;
                if(biotype.equals("three_prime_UTR")) continue;
                if(biotype.equals("five_prime_UTR")) continue;
                if(biotype.equals("CDS")) continue;
                Range e = new Range(chrdb,startdb,enddb);
                rangesList.add(e);
            }
         
            System.out.println("The number of gene ranges is " +  rangesList.size()  );
            
            Rangescp db = new Rangescp(rangesList,"gene");
            db.sortByStartPosition();
            
        for(int i = 0; i<t.getRowNumber(); i++){
        String chrs = t.getCellAsString(i, 0);
        if(chrs.startsWith("bgh")) continue;
        if(chrs.startsWith("chrUn")) continue;
        cnt++;
        chrs = chrs.split("r")[1].split("H")[0];
        int chr = Integer.parseInt(chrs);
        String starts = t.getCellAsString(i, 1); int start = Integer.parseInt(starts);
        String ends = t.getCellAsString(i, 2); int end = Integer.parseInt(ends);
        String count = t.getCellAsString(i, 3).split("_")[1];

        if((db.isInRanges(chr,start)) || (db.isInRanges(chr,end))){
            cntwithingene++;
            bw.write(chrs + "\t" + starts + "\t" + ends + "\t" + count);
            bw.newLine();
        }
        else{
            cntintergenic++;
            bw1.write(chrs + "\t" + starts + "\t" + ends + "\t" + count);
            bw1.newLine();
        }
        }
        System.out.println("The number of reads with chr1H to chr7H is " +  cnt );
        System.out.println("The number of withingene reads  is " +  cntwithingene );
        System.out.println("The number of intergenic reads is " +  cntintergenic );
        br.close();bw.flush();bw.close();bw1.flush();bw1.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public static void main (String[] args){
        new HordeumEntrance();
        System.out.println("\n**********************************" );
        System.out.println("Here is the main class of Hordeum" );
    } 
    
}
