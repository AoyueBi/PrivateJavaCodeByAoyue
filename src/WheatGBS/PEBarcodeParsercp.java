/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGBS;

import analysis.wheat.GBS.*;
import format.table.RowTable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author feilu
 */
public class PEBarcodeParsercp {
    HashMap<String, Set<String>> barcodeTaxonMapR1 = new HashMap<>();
    HashMap<String, Set<String>> barcodeTaxonMapR2 = new HashMap<>();
    String[] barcodesR1 = null;
    String[] barcodesR2 = null;
    String[] taxa = null;
    int maxBarcodeLength = 8;
    public PEBarcodeParsercp (String barcodeFileS) {
        this.buildBarcodeMap(barcodeFileS);
    }
    
    public void buildBarcodeMap (String barcodeFileS) {
        RowTable<String> t = new RowTable(barcodeFileS);
        /*开始处理两端barcode，建立数组barcodesR1和barcodesR2*/
        List<String> barcodeList = t.getColumn(3); /*R1端的barcode列表*/
        Set<String> barcodeSetR1 = new HashSet(barcodeList); /*把R1端的barcode列表去处重复*/
        String[] barcodeArray = barcodeSetR1.toArray(new String[barcodeSetR1.size()]); /*将set集合转化成数组！！！！！<T> T[] toArray(T[] a)*/
        Arrays.sort(barcodeArray); /*将数组排序！！！！Sorts the specified array of objects into ascending order, according
     * to the {@linkplain Comparable natural ordering} of its elements.   我的理解是String[] barcodeeArray是一个对象*/
        barcodesR1 =  barcodeArray; /*把barcodeSetR1转化成String类型数组barcodeArray，然后再进行排序，把排序后的数组赋给barcodesR1*/
        for (int i = 0; i < barcodeArray.length; i++) {
            Set<String> tempSet = new HashSet<String>();
            barcodeTaxonMapR1.put(barcodeArray[i], tempSet); /*把R1端的barcode当作key， tempSet当作value,此时tempSet为空*/
        }
        barcodeList = t.getColumn(4); /*R2端的barcode列表*/
        Set<String> barcodeSetR2 = new HashSet(barcodeList); /*把R2端的barcode列表去处重复*/
        barcodeArray = barcodeSetR2.toArray(new String[barcodeSetR2.size()]);
        Arrays.sort(barcodeArray);
        barcodesR2 = barcodeArray;
        for (int i = 0; i < barcodeArray.length; i++) {
            Set<String> tempSet = new HashSet<String>();
            barcodeTaxonMapR2.put(barcodeArray[i], tempSet); /*把R2端的barcode当作key， tempSet当作value,此时tempSep为空*/
        }
        /*接下来开始处理taxa*/
        taxa = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            String taxon = t.getCell(i, 2);
            taxa[i] = taxon; /*将每行taxon信息赋给数组taxa*/
            String barcode = t.getCell(i, 3); /*将每行1端barcode信息赋给数组barcode*/
            Set<String> tempSet = barcodeTaxonMapR1.get(barcode); /*把barcode作为key,把value赋值给tempSet*/
            tempSet.add(taxon); /*把此时的taxon加入到tempSet中*/
            barcodeTaxonMapR1.put(barcode, tempSet); /*将key和value相关联*/
            barcode = t.getCell(i, 4);
            tempSet = barcodeTaxonMapR2.get(barcode);
            tempSet.add(taxon);
            barcodeTaxonMapR2.put(barcode, tempSet); /*将key和value相关联*/
        }
        Arrays.sort(taxa);
    }
    
    public String getTaxonFromReads (String read1, String read2) {
        int index = Arrays.binarySearch(barcodesR1, read1.substring(0, maxBarcodeLength)); /*在已经排好序的数组barcodeR1中搜索序列前8个碱基是否有匹配的值，
        若有匹配的碱基，即GAATCTA在数组中可以找到，则返回数组的index,此index等于-38，为负*/
        if (index < 0) index = -index-2; /*此时index等于36，说明序列与数组barcodesR1的第37个元素匹配*/
        if (barcodesR1[index].equals(read1.substring(0, barcodesR1[index].length()))) { /*如果我们能从数组中找到与输入序列匹配的碱基，则进行以下操作*/
            Set<String> taxaR1 = barcodeTaxonMapR1.get(barcodesR1[index]); /*在HashMap1中，通过barcode找到对应的taxon，但是一个barcode包含多个taxon名字*/
            if (taxaR1.size() == 1) { /*如果只含有1个taxon的话，将taxaR1从set类型转化为数组类型*/
                return taxaR1.toArray(new String[taxaR1.size()])[0];
            }
            else { /*如果一端barcode含有2个及以上taxon的话*/
                index = Arrays.binarySearch(barcodesR2, read2.substring(0, maxBarcodeLength));
                if (index < 0) index = -index-2;
                if (barcodesR2[index].equals(read2.substring(0, barcodesR2[index].length()))) {
                    Set<String> taxaR2 = barcodeTaxonMapR2.get(barcodesR2[index]);
                    if (taxaR2.size() == 1) {
                        return taxaR2.toArray(new String[taxaR2.size()])[0];
                    }
                    else {
                        taxaR1.retainAll(taxaR2); /*假设taxaR1有2个元素A01 A02 ，taxaR2有3个元素A01 A03 ，则通过retainAll方法保留taxaR1中相同的元素，即A01*/
                        if (taxaR1.size() == 1) {
                            return taxaR1.toArray(new String[taxaR1.size()])[0];
                        }
                        return null;
                    }
                }
            }
        }
        return null;
    }
    
    public String[] getTaxa () {
        String[] anotherTaxa = new String[this.taxa.length];
         /*@param      src      the source array.
     * @param      srcPos   starting position in the source array.
     * @param      dest     the destination array.
     * @param      destPos  starting position in the destination data.
     * @param      length   the number of array elements to be copied.*/
        System.arraycopy(taxa, 0, anotherTaxa, 0, taxa.length);       
        return anotherTaxa;
    }
    
    public HashMap<String, String> getBarcodeTaxaMap () {
        return this.getBarcodeTaxaMap();
    }
}
