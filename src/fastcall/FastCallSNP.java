/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fastcall;

import static cern.jet.math.Arithmetic.factorial;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.set.hash.TByteHashSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.StringReader;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;
import org.apache.commons.math3.stat.inference.ChiSquareTest;


/**
 *
 * @author Fei Lu
 */
public class FastCallSNP {
    int[] chroms = null;
    int[] chromLength = null;
    String[] taxaNames = null;
    String[] bamPaths = null;
    HashMap<String, String[]> taxaBamPathMap = null;
    HashMap<String, String> bamPathPileupPathMap = null;
    ConcurrentHashMap<Integer, Double> factorialMap = new ConcurrentHashMap();
    int maxFactorial = 150;
    Fasta genomeFa = null;
    //A, C, G, T
    byte[] bases = {65, 67, 71, 84};
    //D, I
    byte[] possibleIndel = {68, 73};
    //A, C, D, G, I, T
    byte[] possibleAllele = {65, 67, 68, 71, 73, 84};
    
    
    int binSize = 100000;
    double individualDepthRatioThresh = 0.4; //个体深度率
    double individualThirdAlleleRatioThresh = 0.2; //个体第三个等位基因率
    double segregationPValueThresh = 0.001; //分离P值
    double sequencingErrorRate = 0.1; //测序错误率
    
    
    public FastCallSNP () {
        this.creatFactorialMap();
        this.localTest(1);
    }
    
    public FastCallSNP (String parameterFileS) {
        this.callSNP(parameterFileS);
    }
    
    public void callSNP (String parameterFileS) {
        ArrayList<String> pLineList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(parameterFileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("Author: Fei Lu")) ifOut = true;
            if (!(temp = br.readLine()).equals("Email: flu@genetics.ac.cn; dr.lufei@gmail.com")) ifOut = true;
            if (!(temp = br.readLine()).equals("Homepage: https://plantgeneticslab.weebly.com/")) ifOut = true;
            if (ifOut) {
                System.out.println("Thanks for using FastCall.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String referenceFileS = pLineList.get(0);
        String bamDirS = pLineList.get(1);
        String taxaBamMapFileS = pLineList.get(2);
        int currentChr = Integer.valueOf(pLineList.get(3));
        String vcfDirS = pLineList.get(4);
        
        genomeFa = new Fasta(referenceFileS);
        genomeFa.sortRecordByNameValue(); //对fasta文件进行排序，按照记录的名字。
        chroms = new int[genomeFa.getSeqNumber()];  //对fasta文件进行统计，得出序列的数目，即染色体数目。
        chromLength = new int[genomeFa.getSeqNumber()]; 
        for (int i = 0; i < genomeFa.getSeqNumber(); i++) {
            chroms[i] = Integer.valueOf(genomeFa.getName(i)); //将第一条序列的名字转化为整数类型的值。
            chromLength[i] = genomeFa.getSeqLength(i); //将第一条序列的长度转化为第一号染色体的长度。
        }
        String pileupDirS = new File(new File(vcfDirS).getParent(), "pileup").getAbsolutePath();      
        new File(pileupDirS).mkdir(); //先给出plieup文件的路径，再重新创建文件夹。
        new File(vcfDirS).mkdir(); //创建文件夹vcfDirS
        this.getTaxaBamMap(taxaBamMapFileS); //使用getTaxaBamMap得到一个什么样的结果？taxa的数目，bam文件的数目，每个taxa对应多少个bam文件的HashMap
        File[] bams = new File(bamDirS).listFiles();
        Arrays.sort(bams); //将bam文件进行排序， bams是一个File类型的数组
        this.updateTaxaBamPathMap(bams); //对taxa-bam文件进行更新
        this.creatPileupMap(pileupDirS); //使一个bam文件对应一个pileup文件
        this.creatFactorialMap(); //创建一个i  factorial(i) 的因子图谱
        this.callSNPByChromosome(currentChr, referenceFileS, vcfDirS);
        System.out.println("Variant calling completed");
    }
    
    private void callSNPByChromosome (int currentChr, String referenceFileS, String vcfDirS) {
        int chrIndex = Arrays.binarySearch(chroms, currentChr);
        String chrSeq = genomeFa.getSeq(chrIndex).toUpperCase();
        int regionStart = 1;
        int regionEnd = chrSeq.length();
        this.performPileup(currentChr, regionStart, regionEnd, referenceFileS);
        String outfileS = "chr"+FStringUtils.getNDigitNumber(3, currentChr)+".VCF.txt";
        outfileS = new File (vcfDirS, outfileS).getAbsolutePath(); //VCF输出文件的路径
        int[][] binBound = this.creatBins(currentChr, binSize, regionStart, regionEnd); //建立二维数组的bin， binsize是10万；该条染色体的起始位置信息由上文得知。
        try {
            HashMap<String, BufferedReader> bamPathPileupReaderMap = this.getBamPathPileupReaderMap(); //建立一个HashMap，目的是将一个bam文件对应一个pileup结果的字符缓冲流
            ConcurrentHashMap<BufferedReader, List<String>> readerRemainderMap = this.getReaderRemainderMap(bamPathPileupReaderMap); //通过entryset取出所有value,放入一个空的list中
            BufferedWriter bw = IoUtils.getTextWriter(outfileS); //一个bam文件，写出一个，接下来先写注释信息
            bw.write(this.getAnnotation(referenceFileS)); //注释文件里有换行符
            bw.write(this.getVCFHeader()); //VCF文件的#行
            bw.newLine();
            for (int i = 0; i < binBound.length; i++) {
                long startTimePoint = System.nanoTime();
                int binStart = binBound[i][0];
                int binEnd = binBound[i][1]; //指一个区间，被分成10万bp。是为了将BamPlieup的结果按照每10万个bp来进行一一处理，即10kb，因为？？？？？？？？？？？？？？？？
                ConcurrentHashMap<String, List<List<String>>> bamPileupResultMap = this.getBamPileupResultMap(currentChr, binStart, binEnd, bamPathPileupReaderMap, readerRemainderMap);
                StringBuilder[][] baseSb = this.getPopulateBaseBuilder(binStart, binEnd);
                int[][] depth = this.getPopulatedDepthArray(binStart, binEnd);
                this.fillDepthAndBase(bamPileupResultMap, baseSb, depth, binStart);
                String[][] base = this.getBaseMatrix(baseSb);
                ArrayList<Integer> positionList = this.getPositionList(binStart, binEnd);
                ConcurrentHashMap<Integer, String> posVCFMap = new ConcurrentHashMap((int)((binEnd - binStart + 1)*1.5));
                this.calculateVCF(posVCFMap, positionList, currentChr, binStart, chrSeq, depth, base);
                for (int j = 0; j < positionList.size(); j++) {
                    String vcfStr = posVCFMap.get(positionList.get(j));
                    if (vcfStr == null) continue;
                    bw.write(vcfStr);
                    bw.newLine();
                }
                StringBuilder sb = new StringBuilder();
                sb.append("Bin from ").append(binStart).append(" to ").append(binEnd).append(" is finished. Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
                System.out.println(sb.toString());
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Chromosome " + String.valueOf(currentChr) + " is finished. File written to " + outfileS + "\n");
    }
    
    private void calculateVCF (ConcurrentHashMap<Integer, String> posVCFMap, List<Integer> positionList, int currentChr, int startPos, String chrSeq, int[][] depth, String[][] base) {
        positionList.parallelStream().forEach(position -> {
            int index = position-startPos;
            byte refBase = (byte)(chrSeq.charAt(position-1));
            int baseIndex = Arrays.binarySearch(bases, refBase);
            if (baseIndex < 0) {
                
            }
            else {
                String vcfStr = this.getVCFString(base[index], depth[index], currentChr, position, refBase);
                if (vcfStr != null) {
                    posVCFMap.put(position, this.getVCFString(base[index], depth[index], currentChr, position, refBase));
                }
            }
        });
    }
    
    private String getVCFString (String[] base, int[] depth, int currentChr, int position, byte refBase) {
        TByteArrayList bList;
        boolean ifRecordedDeletion = false;
        TIntHashSet insertionLengthSet = new TIntHashSet();
        TIntHashSet deletionLengthSet = new TIntHashSet();
        int[][] pAlleleCount = new int[base.length][this.possibleAllele.length];
        int[] refDepth = new int[base.length];
        for (int i = 0; i < base.length; i++) {
            bList = new TByteArrayList();
            byte[] ba = base[i].getBytes();
            for (int j = 0; j < ba.length; j++) {
                if (ba[j] == '.') {
                    bList.add(refBase);
                }
                else if (ba[j] == ',') {
                    bList.add(refBase);
                }
                else if (ba[j] == 'A') {
                    bList.add((byte)65);
                }
                else if (ba[j] == 'a') {
                    bList.add((byte)65);
                }
                else if (ba[j] == 'C') {
                    bList.add((byte)67);
                }
                else if (ba[j] == 'c') {
                    bList.add((byte)67);
                }
                else if (ba[j] == 'G') {
                    bList.add((byte)71);
                }
                else if (ba[j] == 'g') {
                    bList.add((byte)71);
                }
                else if (ba[j] == 'T') {
                    bList.add((byte)84);
                }
                else if (ba[j] == 't') {
                    bList.add((byte)84);
                }
                else if (ba[j] == '+') {
                    int endIndex = j+2;
                    for (int k = j+1; k < ba.length; k++) {
                        if (ba[k] > 57) {
                            endIndex = k;
                            break;
                        }
                    }
                    StringBuilder sb = new StringBuilder();
                    for (int k = j+1; k < endIndex; k++) {
                        sb.append((char)ba[k]);
                    }
                    int length = Integer.valueOf(sb.toString());
                    insertionLengthSet.add(length);
                    j+=sb.length();
                    j+=length;
                    if (ba[j-1] == '.' || ba[j-1] == ',') {
                        bList.add((byte)73);
                    }                   
                }
                else if (ba[j] == '-') {
                    int endIndex = j+2;
                    for (int k = j+1; k < ba.length; k++) {
                        if (ba[k] > 57) {
                            endIndex = k;
                            break;
                        }
                    }
                    StringBuilder sb = new StringBuilder();
                    for (int k = j+1; k < endIndex; k++) {
                        sb.append((char)ba[k]);
                    }
                    int length = Integer.valueOf(sb.toString());
                    deletionLengthSet.add(length);
                    j+=sb.length();
                    j+=length;
                    if (ba[j-1] == '.' || ba[j-1] == ',') {
                        bList.add((byte)68);
                    }
                }
                else if (ba[j] == '^') {
                    j++;
                }
                else if (ba[j] == '*') {
                    bList.add(refBase);
                    ifRecordedDeletion = true;
                }
                //N, n, $, >, <
                else {
                    //do nothing
                }
            }
            byte[] taxonBase = bList.toArray();
            for (int j = 0; j < taxonBase.length; j++) {
                int index = Arrays.binarySearch(this.possibleAllele, taxonBase[j]);
                pAlleleCount[i][index]++;
            }
            int altSum = 0;
            for (int j = 0; j < pAlleleCount[i].length; j++) {
                if (this.possibleAllele[j] == refBase) continue;
//                if (this.possibleAllele[j] == 68) continue;
//                if (this.possibleAllele[j] == 73) continue;
                altSum+=pAlleleCount[i][j];
            }
            refDepth[i] = depth[i] - altSum;
        }
        int totalDepth = 0;
        for (int i = 0; i < depth.length; i++) {
            totalDepth+=depth[i];
        }
        
        int[] indelTypeCount = new int[2];
        indelTypeCount[0] = deletionLengthSet.size();
        indelTypeCount[1] = insertionLengthSet.size();
        
        //****************************Filter1 IndelFilter*****************************************
        //Too many indels usually means mis-alignment
        if (indelTypeCount[0]+indelTypeCount[1] > 1) return null;
        //=======================================Filter1=====================================================
        
        //****************************Filter2 Depth_ratio_test*****************************************
        //In any individual, alt allele show up < 2 times, ignore
        //In any individual, alt allele show up 2 times when 2 < depth < 5. Pick up
        //In any individual, alt allele show up > individualDepthRatioThresh, when depth >= 5. Pick up
        //When depth is low, tend to have assembly errors, LTR
        
        TByteHashSet altAlleleSet = new TByteHashSet();
        for (int i = 0; i < pAlleleCount[0].length; i++) {
            if (possibleAllele[i] == refBase) continue;
            for (int j = 0; j < pAlleleCount.length; j++) {
                if (depth[j] < 2) {}
                else if (depth[j] < 5) {
                    if (pAlleleCount[j][i] >=2) {
                        altAlleleSet.add(this.possibleAllele[i]);
                    }
                }
                else {
                    if ((double)pAlleleCount[j][i]/depth[j] > this.individualDepthRatioThresh) {
                        altAlleleSet.add(this.possibleAllele[i]);
                    }
                }
            }
        }
        byte[] altAllele = altAlleleSet.toArray();
        if (altAllele.length == 0) return null;
        Arrays.sort(altAllele);
        //=======================================Filter2=====================================================  
        
        
        int[] altAllele2PAlleleIndex = new int[altAllele.length];
        for (int i = 0; i < this.possibleAllele.length; i++) {
            int index = Arrays.binarySearch(altAllele, this.possibleAllele[i]);
            if (index < 0) continue;
            altAllele2PAlleleIndex[index] = i;
        }   
        int[] altAlleleTotalDepth = new int[altAllele.length];
        for (int i = 0; i < altAllele.length; i++) {
            for (int j = 0; j < base.length; j++) {
                altAlleleTotalDepth[i]+=pAlleleCount[j][altAllele2PAlleleIndex[i]];
            }
        }
        int[] altAlleleDepthDesendingIndex = FArrayUtils.getIndexByDescendingValue(altAlleleTotalDepth);
        int refTotalDepth = 0;
        for (int i = 0; i < refDepth.length; i++) refTotalDepth+=refDepth[i];
        
        //****************************Filter3 third_allele_test************************************************
        //individual should not have the third allele
        if (altAllele.length > 1) {
            for (int i = 0; i < base.length; i++) {
                int[] tempCnt = new int[altAllele.length];
                for (int j = 0; j < altAllele.length; j++) {
                    tempCnt[j] = pAlleleCount[i][altAllele2PAlleleIndex[j]];
                }
                int sum = refDepth[i];
                double[] v = new double[altAllele.length+1];
                for (int j = 0; j < tempCnt.length; j++) v[j] = (double)tempCnt[j]/sum;
                v[v.length-1] = (double)refDepth[i]/sum;
                Arrays.sort(v);
                if (v[v.length-3] > individualThirdAlleleRatioThresh) return null;
            }
        }
        //===========================Filter3=========================================================
        
        //****************************Filter4 Segregation_test*****************************************
        long[] observed = new long[base.length];
        double[] expected = new double[base.length];
        double[] segregationP = new double[altAllele.length];
        ChiSquareTest ct = new ChiSquareTest();
        int cnt = 0;
        for (int i = 0; i < altAllele.length; i++) {
            double r = (double)altAlleleTotalDepth[i]/(refTotalDepth+altAlleleTotalDepth[i]);
            for (int j = 0; j < base.length; j++) {
                observed[j] = pAlleleCount[j][altAllele2PAlleleIndex[i]];
                expected[j] = r;
            }
            segregationP[i] = ct.chiSquareTest(expected, observed);
            if (segregationP[i] > this.segregationPValueThresh) cnt++;
        }
        if (cnt == altAllele.length) return null;
       //===========================Filter4=========================================================
        
        int nonMissingCnt = 0;
        int[] refAndAllelePresence = new int[altAllele.length+1];
        for (int i = 0; i < base.length; i++) {
            if (refDepth[i] != 0) {
                nonMissingCnt++;
                refAndAllelePresence[0]++;
            }
            else {
                for (int j = 0; j < altAllele.length; j++) {
                    if (pAlleleCount[i][altAllele2PAlleleIndex[j]] != 0) {
                        nonMissingCnt++;
                        break;
                    }
                }
            }
            for (int j = 0; j < altAllele.length; j++) {
                if (pAlleleCount[i][altAllele2PAlleleIndex[j]] != 0) {
                   refAndAllelePresence[j+1]++;
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append(currentChr).append("\t").append(position).append("\t.\t").append((char)refBase).append("\t");
        for (int i = 0; i < altAllele.length; i++) sb.append(String.valueOf((char)altAllele[altAlleleDepthDesendingIndex[i]])).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t.\t.\t").append("DP=").append(totalDepth).append(";AD=").append(refTotalDepth);
        for (int i = 0; i < altAllele.length; i++) sb.append(",").append(altAlleleTotalDepth[altAlleleDepthDesendingIndex[i]]);
        sb.append(";NZ=").append(nonMissingCnt).append(";AP=");
        for (int i = 0; i < refAndAllelePresence.length; i++) sb.append(refAndAllelePresence[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append(";PV=");
        for (int i = 0; i < altAllele.length; i++) sb.append(segregationP[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append(";DI=");
        for (int i = 0; i < indelTypeCount.length; i++) sb.append(indelTypeCount[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t").append("GT:AD:PL");
        
        for (int i = 0; i < base.length; i++) {
            int[] dep = new int[altAllele.length+1];
            dep[0] = refDepth[i];
            for (int j = 0; j < altAllele.length; j++) {
                dep[j+1] = pAlleleCount[i][altAllele2PAlleleIndex[altAlleleDepthDesendingIndex[j]]];
            }
            sb.append("\t").append(this.getGenotype(dep));
        }
        return sb.toString();
    }
    
    private ArrayList<Integer> getPositionList (int startPos, int endPos) {
        ArrayList<Integer> positionList = new ArrayList();
        for (int i = startPos; i <= endPos; i++) {
            positionList.add(i);
        }
        return positionList;
    }
    
    private String[][] getBaseMatrix (StringBuilder[][] baseSb) {
        String[][] base = new String[baseSb.length][baseSb[0].length];
        for (int i = 0; i < base.length; i++) {
            for (int j = 0; j < base[0].length; j++) base[i][j] = baseSb[i][j].toString();
        }
        return base;
    }
    
    private void fillDepthAndBase (ConcurrentHashMap<String, List<List<String>>> bamPileupResultMap, StringBuilder[][] baseSb, int[][] depth, int startPos) {
        Set<Map.Entry<String, String[]>> entries = this.taxaBamPathMap.entrySet();
        entries.parallelStream().forEach(e -> {
            String taxa = e.getKey();
            int taxaIndex = Arrays.binarySearch(this.taxaNames, taxa);
            String[] bams = e.getValue();
            
            int count = 0;
            String b = null;
            try {
                for (int i = 0; i < bams.length; i++) {
                    List<List<String>> lines = bamPileupResultMap.get(bams[i]);
                    count = lines.size();
                    b = bams[i];
                    for (int j = 0; j < lines.size(); j++) {
                        List<String> split = lines.get(j);                 
                        if (split.get(2).startsWith("N") || split.get(2).startsWith("n")) continue;
                        int siteIndex = Integer.valueOf(split.get(1)) - startPos;
                        depth[siteIndex][taxaIndex]+=Integer.valueOf(split.get(3));
                        baseSb[siteIndex][taxaIndex].append(split.get(4));
                    }
                }
            }
            catch (Exception ee) {
                System.out.println(b);
                System.out.println(count);
                ee.printStackTrace();
                System.exit(1);
            }
        });
    }
    
    private int[][] getPopulatedDepthArray (int startPos, int endPos) {
        int[][] depth = new int[endPos-startPos+1][this.taxaNames.length];
        return depth;
    }
    
    private StringBuilder[][] getPopulateBaseBuilder (int startPos, int endPos) {
        StringBuilder[][] sbs = new StringBuilder[endPos-startPos+1][this.taxaNames.length];
        for (int i = 0; i < sbs.length; i++) {
            for (int j = 0; j < sbs[0].length; j++) sbs[i][j] = new StringBuilder();
        }
        return sbs;
    }
    
    private ConcurrentHashMap<String, List<List<String>>> getBamPileupResultMap (int currentChr, int binStart, int binEnd, HashMap<String, BufferedReader> bamPathPileupReaderMap, ConcurrentHashMap<BufferedReader, List<String>> readerRemainderMap) {
        ArrayList<String> empty = new ArrayList();
        ConcurrentHashMap<String, List<List<String>>> bamPileupMap = new ConcurrentHashMap(2048);
        List<String> bamList = Arrays.asList(bamPaths);
        bamList.parallelStream().forEach(bamFileS -> {
            ArrayList<List<String>> lineList = new ArrayList();
            BufferedReader br = bamPathPileupReaderMap.get(bamFileS);
            List<String> remainder = readerRemainderMap.get(br);
            boolean flag = false;
            if (remainder.size() == 0) { //一开始，remainder是空的，所以size是0.
                String temp = null;
                try {
                    temp = br.readLine(); //只读第一行，
                }
                catch (Exception e) {}
                if (temp != null) { //对第一行进行拆分，取其位置，和binEnd进行比较，如果大于binEnd，说明不在此区间，将br放入split，否则就将split这一行内容加入lineList
                    List<String> split = FStringUtils.fastSplit(temp, "\t");
                    int currentPos = Integer.valueOf(split.get(1));
                    if (currentPos > binEnd) {
                        readerRemainderMap.put(br, split);
                    }
                    else {
                        lineList.add(split);
                        flag = true;
                    }
                }
            }
            else {
                int currentPos = Integer.valueOf(remainder.get(1)); //
                if (currentPos <= binEnd) {
                    lineList.add(remainder);
                    flag = true;
                    readerRemainderMap.put(br, empty);
                }
            }
            if (flag == true) {
                try {
                    String temp;
                    while ((temp = br.readLine()) != null) {
                        List<String> split = FStringUtils.fastSplit(temp, "\t");
                        int currentPos = Integer.valueOf(split.get(1));
                        if (currentPos < binEnd) {
                            lineList.add(split);
                        }
                        else if (currentPos == binEnd){
                            lineList.add(split);
                            readerRemainderMap.put(br, empty);
                            break;
                        }
                        else {
                            readerRemainderMap.put(br, split);
                            break;
                        }
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            bamPileupMap.put(bamFileS, lineList);
        });
        return bamPileupMap;
    }
    
    private ConcurrentHashMap<BufferedReader, List<String>> getReaderRemainderMap (HashMap<String, BufferedReader> bamPathPileupReaderMap) {
        ArrayList<String> empty = new ArrayList();
        ConcurrentHashMap<BufferedReader, List<String>> readerRemainderMap = new ConcurrentHashMap();
        Set<Map.Entry<String, BufferedReader>> enties = bamPathPileupReaderMap.entrySet(); //通过entrySet()方法将map集合中的映射关系取出（这个关系就是Map.Entry类型）
        enties.stream().forEach(e -> {
            readerRemainderMap.put(e.getValue(), empty);
        });
        return readerRemainderMap;
    }
    
    private HashMap<String, BufferedReader> getBamPathPileupReaderMap () { //目的是将一个bam文件对应一个pileup结果的字符缓冲流
        HashMap<String, BufferedReader> bamPathPileupReaderMap = new HashMap();
        try {
            for (int i = 0; i < this.bamPaths.length; i++) {
                String pileupFileS = this.bamPathPileupPathMap.get(bamPaths[i]);
                File pileupF = new File(pileupFileS);
                if (!pileupF.exists()) { //如果pileup文件不存在，则过滤。
                    BufferedReader br = new BufferedReader(new StringReader(""), 1024);
                    bamPathPileupReaderMap.put(bamPaths[i], br);
                    continue;
                }
                BufferedReader br = new BufferedReader (new FileReader(pileupF), 1024);
                bamPathPileupReaderMap.put(bamPaths[i], br);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bamPathPileupReaderMap;
    }
    
    private void performPileup (int currentChr, int startPos, int endPos, String referenceFileS) {
        System.out.println("Pileup is being performed on chromosome "+String.valueOf(currentChr)+" from "+String.valueOf(startPos)+" to "+String.valueOf(endPos));
        long timeStart = System.nanoTime();
        List<String> bamList = Arrays.asList(bamPaths);
        LongAdder counter = new LongAdder(); //如果程序内有高度的竞争，大量的线程访问同一个原子值，可以使用 LongAdder 和 LongAccumulator，这个类是 Java 8 提供用于在高度竞争环境下替代 AtomicLong 的。
        bamList.parallelStream().forEach(bamFileS -> {
            String pileupFileS = this.bamPathPileupPathMap.get(bamFileS);
            StringBuilder sb = new StringBuilder();
            sb.append("samtools mpileup -A -B -q 30 -Q 10 -f ").append(referenceFileS).append(" ").append(bamFileS).append(" -r "); //-r区域
            //-A能够帮助我们不跳过哪些异常的read pairs; -B帮助我们降低由错配引起的假阳性SNPs; -q 限定最小比对质量为30; -Q限定最小碱基质量为10;-f为参考基因组的fa格式或bgzip压缩过且带有fai文件
            sb.append(currentChr).append(":").append(startPos).append("-").append(endPos).append(" -o ").append(pileupFileS);
            String command = sb.toString();
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            counter.increment();
            int cnt = counter.intValue();
            if (cnt%10 == 0) System.out.println("Pileuped " + String.valueOf(cnt) + " bam files. Total: " + String.valueOf(this.bamPaths.length));
        });
        System.out.println("Pileup is finished. Time took " + Benchmark.getTimeSpanMinutes(timeStart) + " mins");
    }
     
    private int[][] creatBins (int currentChr, int binSize, int regionStart, int regionEnd) {
        int[][] binBound = FArrayUtils.getSubsetsIndicesBySubsetSize(regionEnd-regionStart+1, binSize);
        for (int i = 0; i < binBound.length; i++) {
            binBound[i][0]+=regionStart;
            binBound[i][1]+=regionStart;
            binBound[i][1]--;
        }
        System.out.println("SNP calling will performed on chromosome "+String.valueOf(currentChr)+" from "+String.valueOf(regionStart)+" to "+String.valueOf(regionEnd));
        System.out.println("Chromosome " + String.valueOf(currentChr) +"  is devided into " + String.valueOf(binBound.length) + " bins. Bin size: " + String.valueOf(this.binSize) + " bp");
        return binBound;
    }
    
    private String getVCFHeader () {
        StringBuilder sb = new StringBuilder();
        sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (int i = 0; i < taxaNames.length; i++) {
            sb.append("\t").append(taxaNames[i]);
        }
        return sb.toString();
    }
    
    
    private void creatFactorialMap () {
        for (int i = 0; i < this.maxFactorial+1; i++) {
            this.factorialMap.put(i, factorial(i));
        }
    }
    
    public String getGenotype (int[] cnt) {
        int n = cnt.length*(cnt.length+1)/2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum+=cnt[i];
        if (sum == 0) return "./.";
        else if (sum > this.maxFactorial) {
            double portion = (double)this.maxFactorial/sum;
            for (int i = 0; i < cnt.length; i++) {
                cnt[i] = (int)(cnt[i]*portion);
            }
            sum = this.maxFactorial;
        }
        double coe = this.factorialMap.get(sum);
        for (int i = 0; i < cnt.length; i++) coe = coe/this.factorialMap.get(cnt[i]);
        double max = Double.MAX_VALUE;
        int a1 = 0;
        int a2 = 0;
        for (int i = 0; i < cnt.length; i++) {
            for (int j = i; j < cnt.length; j++) {
                int index = (j*(j+1)/2)+i;
                double value = Double.MAX_VALUE;
                if (i == j) {
                    value = -Math.log10(coe*Math.pow((1-0.75*this.sequencingErrorRate), cnt[i])*Math.pow(this.sequencingErrorRate/4, (sum-cnt[i])));
                }
                else {
                    value = -Math.log10(coe*Math.pow((0.5-this.sequencingErrorRate/4), cnt[i]+cnt[j])*Math.pow(this.sequencingErrorRate/4, (sum-cnt[i]-cnt[j])));
                }
                if (value < max) {
                    max = value;
                    a1 = i;
                    a2 = j;
                }
                likelihood[index] = (int)Math.round(value);
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append(a1).append("/").append(a2).append(":");
        for (int i = 0; i < cnt.length; i++) sb.append(cnt[i]).append(",");
        sb.deleteCharAt(sb.length()-1); sb.append(":");
        for (int i = 0; i < likelihood.length; i++) sb.append(likelihood[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
    
    public int[] getGTLikelihood (int[] cnt) {
        int n = cnt.length*(cnt.length+1)/2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum+=cnt[i];
        double coe = factorial(sum);
        for (int i = 0; i < cnt.length; i++) coe = coe/factorial(cnt[i]);
        for (int i = 0; i < cnt.length; i++) {
            for (int j = i; j < cnt.length; j++) {
                int index = (j*(j+1)/2)+i;
                if (i == j) {
                    likelihood[index] = (int)Math.round(-Math.log10(coe*Math.pow((1-0.75*this.sequencingErrorRate), cnt[i])*Math.pow(this.sequencingErrorRate/4, (sum-cnt[i]))));
                }
                else {
                    likelihood[index] = (int)Math.round(-Math.log10(coe*Math.pow((0.5-this.sequencingErrorRate/4), cnt[i]+cnt[j])*Math.pow(this.sequencingErrorRate/4, (sum-cnt[i]-cnt[j]))));
                }
            }
        }
        return likelihood;
    }
    
    private void creatPileupMap (String pileupDirS) { //一个bam文件，对应一个pileup文件
        bamPathPileupPathMap = new HashMap();
        Set<Entry<String, String[]>> entries = taxaBamPathMap.entrySet();
        for (Entry<String, String[]> e : entries) {
            String[] bams = e.getValue();
            for (String bam : bams) {
                String pileupFileS = new File (pileupDirS, new File(bam).getName().replaceFirst(".bam", ".pileup.txt")).getAbsolutePath();
                bamPathPileupPathMap.put(bam, pileupFileS);
            }
        }
    }
      
    private void updateTaxaBamPathMap (File[] bams) {
        String bamDirS = bams[0].getParent(); //存放真实bam文件的路径
        String[] existingBam = new String[bams.length]; //真实文件中bam的个数
        for (int i = 0; i < bams.length; i++) existingBam[i] = bams[i].getName(); //遍历真实bam文件，返回的是真实bam文件的名字
        Arrays.sort(existingBam); //对真实bam文件进行排序
        HashSet<String> existingTaxaSet = new HashSet(); //********先建立一个existingTaxaSet集合，这个集合是空的，以后用来存放现存的taxa？？？
        HashMap<String, String[]> updatedTaxaBamMap = new HashMap(); //**建立更新后的updatedTaxaBamMap HashMap，用于一个taxa对应多个bam文件
        int cnt = 0;
        ArrayList<String> pathList = new ArrayList();  //*********建立一个pathList，用来存放真实bam文件的全路径集合。
        for (int i = 0; i < taxaNames.length; i++) {   //对一直文件列表中的taxa名字进行循环，根据taxa名字，找到对应的bam列表。
            String[] bamNames = taxaBamPathMap.get(taxaNames[i]); //bamNames用来存放taxa对应的原本计划中的bam文件集合
            ArrayList<String> bamPathList = new ArrayList(); //*********建立一个bamPathList集合，用来存放？？？？？？？
            for (int j = 0; j < bamNames.length; j++) { //在每个taxa下循环，找到bam集合后，再进行bam的循环。
                int index = Arrays.binarySearch(existingBam, bamNames[j]);  //库为真实bam文件的个数，bamNames[j]为参数列表里的bam文件。
                if (index < 0) continue; //去掉列表里存在，但在真实库里不存在的bam；
                String path = new File(bamDirS,bamNames[j]).getAbsolutePath(); //真实bam的路径
                bamPathList.add(path); //存放真实bam的路径名
                pathList.add(path); //存放真实bam的路径名
                existingTaxaSet.add(taxaNames[i]); //存放taxa的集合，这里的意思是，必须有真实bam文件存在，这个taxa才是真实存在的，若真实bam文件没有，则列表里的taxa也不会加入到集合中去。
            }
            if (bamPathList.isEmpty()) continue; //如果taxa对于的bam文件为空，则略去这个taxa
            bamNames = bamPathList.toArray(new String[bamPathList.size()]); //此时bamNames存放的是真实bam文件的带路径的文件名
            Arrays.sort(bamNames);
            updatedTaxaBamMap.put(taxaNames[i], bamNames); //此时的HashMap对应的是 一个taxa,一群带bam路径名的数组
            cnt+=bamNames.length;
        }
        String[] updatedTaxaNames = existingTaxaSet.toArray(new String[existingTaxaSet.size()]); //有bam文件的taxa集合
        Arrays.sort(updatedTaxaNames);
        taxaNames = updatedTaxaNames;
        taxaBamPathMap = updatedTaxaBamMap;
        this.bamPaths = pathList.toArray(new String[pathList.size()]);
        Arrays.sort(bamPaths);
        System.out.println("Actual taxa number:\t"+String.valueOf(taxaNames.length));
        System.out.println("Actual bam file number:\t"+String.valueOf(cnt));
        System.out.println();
    }
    
    private void getTaxaBamMap (String taxaBamMapFileS) {
        this.taxaBamPathMap = new HashMap(); // HashMap<String, String[]> taxaBamPathMap 先为空，后期再进行创建！new一个新的对象
        try {
            BufferedReader br = IoUtils.getTextReader(taxaBamMapFileS);
            String temp;
            ArrayList<String> taxaList = new ArrayList(); //所以taxa名字的集合
            int nBam = 0;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                taxaList.add(tem[0]);
                String[] bams = new String[tem.length-1] ; //一个taxa对应多个bam文件，组成一个数组
                for (int i = 0; i < bams.length; i++) bams[i] = tem[i+1];
                Arrays.sort(bams); //bam数组进行排序
                taxaBamPathMap.put(tem[0], bams); //每一个taxa对应一组bam文件，组成一个HashMap
                nBam+=bams.length; //统计所有bam文件的个数
            }
            taxaNames = taxaList.toArray(new String[taxaList.size()]);  //将taxa结合转化为数组
            Arrays.sort(taxaNames);
            System.out.println("Created TaxaBamMap from" + taxaBamMapFileS);
            System.out.println("Taxa number:\t"+String.valueOf(taxaNames.length));
            System.out.println("Bam file number in TaxaBamMap:\t"+String.valueOf(nBam));
            System.out.println();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    //********************************Local test section***********************************   
    public void localTest (int currentChr) {
        String taxaBamMapFileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\GenotypingPipelineTest\\source\\taxaBamMap.txt";
        String pileupDirS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\GenotypingPipelineTest\\pileup\\";
        String vcfDirS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\GenotypingPipelineTest\\vcfFast\\";
        String bamList = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\GenotypingPipelineTest\\bam.list.txt" ;
        String referenceFileS = "M:\\Database\\cassavaReference\\genome\\Manihot esculenta\\cassavaV6_chr01.fa";
        int startPos = 1;
        int endPos = 100000;
        int binSize = 10000;
        Fasta f = new Fasta(referenceFileS);
        String chrSeq = f.getSeq(0).toUpperCase();
        this.getTaxaBamMap(taxaBamMapFileS);
        this.updateTaxaBamPathMap(bamList);
        this.creatPileupMap(pileupDirS);
        this.callSNPByChromosomeTest(currentChr, vcfDirS, referenceFileS, chrSeq, startPos, endPos, binSize);
    }
    
    private void callSNPByChromosomeTest (int currentChr, String vcfDirS, String referenceFileS, String chrSeq, int startPos, int endPos, int binSize) {
        int regionStart = startPos;
        int regionEnd = endPos;
        
        String outfileS = "chr"+FStringUtils.getNDigitNumber(3, currentChr)+".VCF.txt";
        outfileS = new File (vcfDirS, outfileS).getAbsolutePath();
        int[][] binBound = this.creatBins(currentChr, binSize, regionStart, regionEnd);
        try {
            HashMap<String, BufferedReader> bamPathPileupReaderMap = this.getBamPathPileupReaderMap();
            ConcurrentHashMap<BufferedReader, List<String>> readerRemainderMap = this.getReaderRemainderMap(bamPathPileupReaderMap);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(this.getAnnotation(referenceFileS));
            bw.write(this.getVCFHeader());
            bw.newLine();
            for (int i = 0; i < binBound.length; i++) {
                long startTimePoint = System.nanoTime();
                int binStart = binBound[i][0];
                int binEnd = binBound[i][1];
                ConcurrentHashMap<String, List<List<String>>> bamPileupResultMap = this.getBamPileupResultMap(currentChr, binStart, binEnd, bamPathPileupReaderMap, readerRemainderMap);
                StringBuilder[][] baseSb = this.getPopulateBaseBuilder(binStart, binEnd);
                int[][] depth = this.getPopulatedDepthArray(binStart, binEnd);
                this.fillDepthAndBase(bamPileupResultMap, baseSb, depth, binStart);
                String[][] base = this.getBaseMatrix(baseSb);
                ArrayList<Integer> positionList = this.getPositionList(binStart, binEnd);
                ConcurrentHashMap<Integer, String> posVCFMap = new ConcurrentHashMap((int)((binEnd - binStart + 1)*1.5));
                this.calculateVCF(posVCFMap, positionList, currentChr, binStart, chrSeq, depth, base);
                for (int j = 0; j < positionList.size(); j++) {
                    String vcfStr = posVCFMap.get(positionList.get(j));
                    if (vcfStr == null) continue;
                    bw.write(vcfStr);
                    bw.newLine();
                }
                StringBuilder sb = new StringBuilder();
                sb.append("Bin from ").append(binStart).append(" to ").append(binEnd).append(" is finished. Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
                System.out.println(sb.toString());
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Chromosome " + String.valueOf(currentChr) + " is finished. File written to " + outfileS + "\n");
    }
    
    private String getAnnotation (String referenceFileS) {
        StringBuilder sb = new StringBuilder();
        sb.append("##fileformat=VCFv4.1\n");
        SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
        Date dt = new Date();
        String S = sdf.format(dt);
        sb.append("##fileDate=").append(S.split(" ")[0]).append("\n");
        sb.append("##reference=").append(referenceFileS).append("\n");
        sb.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"").append("Total depth").append("\">\n");
        sb.append("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"").append("Total allelelic depths in order listed starting with REF").append("\">\n");
        sb.append("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"").append("Number of individuals with alleles present").append("\">\n");
        sb.append("##INFO=<ID=AP,Number=2+,Type=Integer,Description=\"").append("Number of individuals in which an allele is present").append("\">\n");
        sb.append("##INFO=<ID=PV,Number=1+,Type=Float,Description=\"").append("Segreagation test P-Value of alternative alleles").append("\">\n");
        sb.append("##INFO=<ID=DI,Number=2,Type=Integer,Description=\"").append("Number of deletion and insertion type").append("\">\n");
        sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"").append("Genotype").append("\">\n");
        sb.append("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"").append("Allelic depths for the reference and alternate alleles in the order listed").append("\">\n");
        sb.append("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"").append("Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles").append("\">\n");
        return sb.toString();
    }
    
    private void updateTaxaBamPathMap (String bamListFileS) {
        Table t = new Table (bamListFileS);
        String[] existingBam = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) existingBam[i] = t.content[i][0];
        Arrays.sort(existingBam);
        HashSet<String> existingTaxaSet = new HashSet();
        HashMap<String, String[]> updatedTaxaBamMap = new HashMap();
        ArrayList<String> pathList = new ArrayList();
        int cnt = 0;
        for (int i = 0; i < taxaNames.length; i++) {
            String[] bamNames = taxaBamPathMap.get(taxaNames[i]);
            ArrayList<String> bamList = new ArrayList();
            for (int j = 0; j < bamNames.length; j++) {
                int index = Arrays.binarySearch(existingBam, bamNames[j]);
                if (index < 0) continue;
                bamList.add(bamNames[j]);
                existingTaxaSet.add(taxaNames[i]);
                pathList.add(bamNames[j]);
            }
            if (bamList.isEmpty()) continue;
            bamNames = bamList.toArray(new String[bamList.size()]);
            Arrays.sort(bamNames);
            updatedTaxaBamMap.put(taxaNames[i], bamNames);
            cnt+=bamNames.length;
        }
        String[] updatedTaxaNames = existingTaxaSet.toArray(new String[existingTaxaSet.size()]);
        Arrays.sort(updatedTaxaNames);
        taxaNames = updatedTaxaNames;
        taxaBamPathMap = updatedTaxaBamMap;
        this.bamPaths = pathList.toArray(new String[pathList.size()]);
        Arrays.sort(bamPaths);
        System.out.println("Actual taxa number:\t"+String.valueOf(taxaNames.length));
        System.out.println("Actual bam file number:\t"+String.valueOf(cnt));
        System.out.println();
    }  
//************************************************************************************************************  
    
    public static void main (String[] args) {
        System.out.println("This is the entrance of FastCall");
        new FastCallSNP (args[0]);
        //new FastCallSNP (); //for test
    }
}
