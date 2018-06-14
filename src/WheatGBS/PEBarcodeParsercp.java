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
        List<String> barcodeList = t.getColumn(3);
        Set<String> barcodeSetR1 = new HashSet(barcodeList);
        String[] barcodeArray = barcodeSetR1.toArray(new String[barcodeSetR1.size()]);
        Arrays.sort(barcodeArray);
        barcodesR1 =  barcodeArray; 
        for (int i = 0; i < barcodeArray.length; i++) {
            Set<String> tempSet = new HashSet<String>();
            barcodeTaxonMapR1.put(barcodeArray[i], tempSet);
        }
        barcodeList = t.getColumn(4);
        Set<String> barcodeSetR2 = new HashSet(barcodeList);
        barcodeArray = barcodeSetR2.toArray(new String[barcodeSetR2.size()]);
        Arrays.sort(barcodeArray);
        barcodesR2 = barcodeArray;
        for (int i = 0; i < barcodeArray.length; i++) {
            Set<String> tempSet = new HashSet<String>();
            barcodeTaxonMapR2.put(barcodeArray[i], tempSet);
        }
        taxa = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            String taxon = t.getCell(i, 2);
            taxa[i] = taxon;
            String barcode = t.getCell(i, 3);
            Set<String> tempSet = barcodeTaxonMapR1.get(barcode);
            tempSet.add(taxon);
            barcodeTaxonMapR1.put(barcode, tempSet);
            barcode = t.getCell(i, 4);
            tempSet = barcodeTaxonMapR2.get(barcode);
            tempSet.add(taxon);
            barcodeTaxonMapR2.put(barcode, tempSet);
        }
        Arrays.sort(taxa);
    }
    
    public String getTaxonFromReads (String read1, String read2) {
        int index = Arrays.binarySearch(barcodesR1, read1.substring(0, maxBarcodeLength));
        if (index < 0) index = -index-2;
        if (barcodesR1[index].equals(read1.substring(0, barcodesR1[index].length()))) {
            Set<String> taxaR1 = barcodeTaxonMapR1.get(barcodesR1[index]);
            if (taxaR1.size() == 1) {
                return taxaR1.toArray(new String[taxaR1.size()])[0];
            }
            else {
                index = Arrays.binarySearch(barcodesR2, read2.substring(0, maxBarcodeLength));
                if (index < 0) index = -index-2;
                if (barcodesR2[index].equals(read2.substring(0, barcodesR2[index].length()))) {
                    Set<String> taxaR2 = barcodeTaxonMapR2.get(barcodesR2[index]);
                    if (taxaR2.size() == 1) {
                        return taxaR2.toArray(new String[taxaR2.size()])[0];
                    }
                    else {
                        taxaR1.retainAll(taxaR2);
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
        System.arraycopy(taxa, 0, anotherTaxa, 0, taxa.length);
        return anotherTaxa;
    }
    
    public HashMap<String, String> getBarcodeTaxaMap () {
        return this.getBarcodeTaxaMap();
    }
}
