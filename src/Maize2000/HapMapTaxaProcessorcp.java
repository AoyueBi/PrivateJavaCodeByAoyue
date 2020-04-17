/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Maize2000;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;


/**
 *
 * @author feilu
 */
public class HapMapTaxaProcessorcp {
    
    public HapMapTaxaProcessorcp () {
        //this.mkSampleTaxaMap();
        String infileS1 = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/taxaNameMap.txt";
        /*t.txt 文件格式是 Taxa	BamPath
K13HGL0039	/public-supool/home/biaoyue/output_data_rmdup_350bplibrary/K13HGL0039.sorted.rmdup.bam
K13HGL0105	/public-supool/home/biaoyue/output_data_rmdup_350bplibrary/K13HGL0105.sorted.rmdup.bam
        
        taxaBam.txt文件格式是 Taxa	BamPath
Shuang741	/public-supool/home/biaoyue/output_data_rmdup_350bplibrary/K13HGL0039.sorted.rmdup.bam
Q1261	/public-supool/home/biaoyue/output_data_rmdup_350bplibrary/K13HGL0105.sorted.rmdup.bam
*/
        String infileS2 = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/t.txt";
        String outfileS = "/Users/feilu/Desktop/taxaBam.txt";
        //updateTaxaBamMap(infileS1, infileS2, outfileS);
       // this.testTaxaDuplicates();
        this.iftest();
    }
    
    public void iftest(){
        int[] cnt = {2,0,1,1};
        for(int i = 0; i< 4; i++){
            if (cnt[i] < 2) continue;
            System.out.println(cnt[i]);          
        }
    }
    
    private void mkSampleTaxaMap () {
        /*initialTaxaNameMap.txt是 测序样本编号	材料名称    K16HL0359	P39/su	 */
        /*taxaNameMap.txt是 --	-- K16HL0359	P39/su  P39_su	*/
        String infileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/initialTaxaNameMap.txt";
        String outfileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/taxaNameMap.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = new ArrayList(); 
        for (int i = 0; i < t.getRowNumber(); i++) {
            /*currentName是指拿到的taxa名字*/
            String currentName = t.getCell(i, 1);
            if (!isStandardTaxonName(currentName)) {
                currentName = this.getStandardTaxonName(currentName);
            }
            l.add(currentName);
        }
        t.insertColumn("Taxa", 2, l);
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
    
    private String getStandardTaxonName (String taxonName) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < taxonName.length(); i++) {
            int c = (int)taxonName.charAt(i);
            if ((47 < c && c < 58) || ((64 < c && c < 91)) || (96 < c && c < 123) || c == 45 || c ==95) {
                sb.append(taxonName.charAt(i));
            }
            else {
                if (sb.length() != 0) {
                    sb.append("_"); /*把不符合标准的字符，比如空格和斜杠，转化为下标横*/
                }
            }
        }
        return sb.toString(); /*最后返回一个标准的名字*/
    }
    
    private boolean isStandardTaxonName (String taxonName) {
        int cnt = 0;
        for (int i = 0; i < taxonName.length(); i++) {
            /*将一个taxaName的每个字符都拿来判断，charAt() 方法用于返回指定索引处的字符。索引范围为从 0 到 length() - 1。*/
            int c = (int)taxonName.charAt(i);
            if ((47 < c && c < 58) || ((64 < c && c < 91)) || (96 < c && c < 123) || c == 45 || c ==95) {
                cnt++; /*表明taxa，名字的第一个字符判断完毕，是符合标准的就加1，不符合标准的，就返回false，说明命名方式有问题，需要重新改正*/
            }
        }
        if (cnt == taxonName.length()) return true;
        return false;
    }
    
    public static void updateTaxaBamMap (String infileS1, String infileS2, String outfileS) {
        /*思路：将ID - Taxon读进表格，并建立HashMap，将 ID - filepath 也读进表格，并用setCell方法将此表格中的第一列ID替换成Taxa信息，
        重新写进一个新的表格中去，写进新表格的方法为writeTextTable*/
        RowTable<String> t = new RowTable<>(infileS1);
        HashMap<String, String> sampleTaxaMap = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            sampleTaxaMap.put(t.getCell(i, 0), t.getCell(i, 2));
        }
        t = new RowTable<>(infileS2);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String c = t.getCell(i, 0);
            t.setCell(i, 0, sampleTaxaMap.get(c));
        }
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
}
