/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatUnclass;

import format.table.RowTable;
import format.tree.Newick;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class TreeFromShusong {
    public TreeFromShusong(){
        //this.addCol();
        //this.convertrow_Vol();
        //this.selectedTaxaFromNewickTree();
        this.convertTaxaNameLutoZheng();
        //this.findNA();
        //this.findNANA();
        //this.colorSelectedTaxaOnTree();
        //this.addseries();
    }
    
    public void addseries(){
        for(int i = 1; i < 21; i++){
            System.out.println(i);
            System.out.println(i);
            System.out.println(i);
            
        }
    }
    
    public void colorSelectedTaxaOnTree () {
        String infileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/source/Tree_chr10000sites_432taxa.xml";
        String outfileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/Tree_chr10000sites_432taxa_colorSel.xml";
        String taxaFileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/selectedTaxa.txt";
        String addS = "<property ref=\"style:font_color\" datatype=\"xsd:token\" applies_to=\"node\">#ff0000</property>";
        RowTable<String> t = new RowTable(taxaFileS);
        List<String> taxaList = t.getColumn(0);
        Collections.sort(taxaList);
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.contains("<name>")) {
                    String currentName = temp.split(">")[1];
                    currentName = currentName.split("<")[0];
                    bw.write(temp);
                    bw.newLine();
                    bw.write(br.readLine());
                    bw.newLine();
                    if (Collections.binarySearch(taxaList, currentName) < 0) {
                        
                    }
                    else {
                        bw.write(addS);
                        bw.newLine();
                    } 
                }
                else {
                    bw.write(temp);
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
    }
    
    /**
     * 从Taxa381中找到不在344材料里的taxa 40个
     */
    public void findNANA(){
        String infileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/001DataBase/luTaxa344name_convertZheng.txt";
        String dbFileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/001DataBase/zhengTaxa381name.txt";
        String outfileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/001DataBase/NA40name.txt";
        List<String> l = new ArrayList();
        String temp = null;
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            //String temp = null;
            while((temp = br.readLine()) != null){
                String name = PStringUtils.fastSplit(temp).get(0);
                l.add(name);
            }
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.out.println(temp);
            System.exit(1);
        }
        String[] taxa344 = l.toArray(new String[l.size()]);
        Arrays.sort(taxa344);
        
        try{
            BufferedReader br = IOUtils.getTextReader(dbFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp2 = br.readLine();
            while((temp2 = br.readLine()) != null){
                int index = Arrays.binarySearch(taxa344, temp2);
                if(index >= 0) continue;
                bw.write(temp2);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * 找到344份材料中，有3个样不在381里面。
     */
    public void findNA(){
        String infileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/001DataBase/luTaxa344name_convertZheng.txt";
        String dbFileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/001DataBase/zhengTaxa381name.txt";
        String outfileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/001DataBase/NAname.txt";
        RowTable<String> t = new RowTable<>(dbFileS);
        List<String> taxa381List = t.getColumn(0);
        String[] taxa381 = taxa381List.toArray(new String[taxa381List.size()]);
        Arrays.sort(taxa381);
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while((temp = br.readLine()) != null){
                String queryName = PStringUtils.fastSplit(temp).get(0);
                int index = Arrays.binarySearch(taxa381, queryName);
                if(index >= 0) continue;
                bw.write(queryName);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void convertTaxaNameLutoZheng(){
        String infileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/selected60Taxa_from70.txt";
        String dbFileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/001DataBase/DataBase.txt";
        String outfileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/selected60Taxa_from70_conZheng.txt";
        RowTable<String> t = new RowTable<>(dbFileS);
        HashMap<String, String> hm = new HashMap<>();
        for(int i = 0; i<t.getRowNumber(); i++){
            hm.put(t.getCellAsString(i, 1), t.getCellAsString(i, 0));
        }
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            bw.write(temp);bw.newLine();
            while((temp = br.readLine()) != null){
                String key = PStringUtils.fastSplit(temp).get(0);
                bw.write(hm.get(key)); 
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
    public void selectedTaxaFromNewickTree () {
        String infileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/001DataBase/luTaxa344-3.newick";
        String outfileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/selected70Taxa.txt";
        //String[] vip = {"Jing2416", "mPH6WC", "Huangzaosi", "B73", "Mo17Ht"};     
        int n = 70;
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            br.close();
            //temp = "((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154)";
            Newick nwk = new Newick (temp);
            List<String> taxaList = nwk.selectTaxaWithMaxDiversity(n);
            Set<String> taxaSet = new HashSet(taxaList);
//            for (int i = 0; i < vip.length; i++) {
//                taxaSet.add(vip[i]);
//            }
            taxaList = new ArrayList(taxaSet);
            Collections.sort(taxaList);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa");
            bw.newLine();
            for (int i = 0; i < taxaList.size(); i++) {
                bw.write(taxaList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void convertrow_Vol(){
        String infileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/Taxa381_row.txt";
        String outfileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/Taxa381_vol.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while((temp = br.readLine()) != null){
                List<String> l = PStringUtils.fastSplit(temp," ");
                for(int i=0; i < l.size(); i++){
                    bw.write(l.get(i));
                    bw.newLine();
                }
               
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    public void addCol(){
        String infileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/selected60Taxa_from70.txt";
        String outfileS = "/Users/Aoyue/Documents/Tree_ShusongZheng/selected60Taxa_from70_col.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            while((temp = br.readLine()) != null){
                
                bw.write(temp);bw.write("\tred");
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
    public static void main (String[] args){
        new TreeFromShusong();
        
    }  
    
}
