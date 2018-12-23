/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatBWA;

import WheatGBS.DataProcessor;
import WheatGBS.PlateAndID;
import format.table.RowTable;
import format.table.TableInterface;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class WheatBWA {
    public WheatBWA(){
        //this.mkMd5();
        //this.checkMd5();
        //this.fastQC();
        //this.alignBWA();
        //this.listSpecificalFiles();
        //this.testspilt();
        //this.samtoolsSort();
        //this.samtoolsMerge();
        this.pipeline();
        
    }
    
    /**
     * 本程序适合于一个样品的深测序质控；每个脚本单独运行完毕后，开始运行下一个脚本。
     * 输出脚本目录：ScriptParentS
     * 数据输出父目录：HPCS
     * Name列表：infileS
     * 原始数据文件夹：rawdataDirS
     */
    
    public void pipeline(){
        
        String parentFileS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/a-output/"; //在本地主机上 建立工作文件夹
        String HPCS = "/data2/aoyue/a-output/"; //集群工作父目录
        String ScriptParentS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/script/"; //脚本生成文件 父目录
        String bamS = "bamfile/";  //子文件夹
        String mergeS = "mergefile/"; //子文件夹
        String statS = "statsfile/"; //子文件夹
        String testS = "testbamfile/"; //子文件夹
        
        /*********************** 建立工作区(在本地路径中) **************************/
        new File(parentFileS).mkdirs();
        new File(parentFileS,bamS).mkdirs();
        new File(parentFileS,mergeS).mkdirs();
        new File(parentFileS,statS).mkdirs();
        new File(parentFileS,testS).mkdirs();
        
        /*********************** 建立工作区(在HPC中) **************************/
        String mkdirScriptS = ScriptParentS +"000_mkdir.sh";
        try{
            BufferedWriter bw = IOUtils.getTextWriter(mkdirScriptS);
            bw.write("mkdir " + HPCS);bw.newLine();
            bw.write("mkdir " + HPCS + bamS);bw.newLine();
            bw.write("mkdir " + HPCS + mergeS);bw.newLine();
            bw.write("mkdir " + HPCS + statS);bw.newLine();
            bw.write("mkdir " + HPCS + testS);bw.newLine();
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        String infileS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/001_script/SM8.t.txt"; // 建立Name列表
        //String infileS = "/Users/Aoyue/Documents/test2M.txt"; //测试数据
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        Collections.sort(namelist);
        String[] names = namelist.toArray(new String[namelist.size()]);
        Arrays.sort(names);
        
        
        String bwaThread = "8";  //根据实际情况修改，使用的线程数
        String indexFileS = "/data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz"; //比对参考基因组index library
        //String rawdataDirS = "/data2/aoyue/a-wheatRawdata1219";
        String rawdataDirS = "/data2/aoyue/test/"; //原始数据文件夹
        String bamfileDirS = HPCS + bamS; 
        String bwaScriptS = ScriptParentS +"001_bwa.sh";
        
        
        String sortScriptS = ScriptParentS + "002_sort.sh";
        String sortMemory = "10G"; //根据实际情况修改，每个线程使用的最大内存
        String sortThread = "10"; //根据实际情况修改，使用的线程数
        
        String mergefileDirS = HPCS + mergeS;
        String mergeName = "mergeWheat" + String.valueOf(names.length) + "SM.bam";
        String mergeScriptS = ScriptParentS + "003_merge.sh";
        
        String RefSeqS = "/data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa"; //统计stat需要解压的ref fa文件
        String bamQCScriptS = ScriptParentS + "004_bamQC.sh";
        String statfileDirS = HPCS + statS;
        
        
        
        
        
        /**************************** bwa **************************************/
        try {
            BufferedWriter bw = IOUtils.getTextWriter(bwaScriptS);
            for (int i = 0; i < names.length; i++) {
            //bwa mem -t 100 -R '@RG\tID:HRV-L1\tPL:illumina\tSM:HRV-L1\tLB:HRV-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD/HRV-L1_1.fq.gz /data2/sharedData/Jiao/ABD/HRV-L1_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HM.pe.bam && echo "** bwa mapping done **" &
                StringBuilder sb = new StringBuilder("bwa mem -t ");
                sb.append(bwaThread).append(" -R ").append("'@RG\\tID:").append(names[i]).append("\\t").append("PL:illumina").append("\\t").append("SM:").append(names[i]).append("\\t").append("LB:").append(names[i]).append("' ");
                sb.append(indexFileS).append(" ");
                sb.append(new File(rawdataDirS, names[i]+"_1.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(rawdataDirS, names[i]+"_2.fq.gz").getAbsolutePath());
                sb.append(" | ").append("samtools view -S -b - > ").append(new File(bamfileDirS, names[i]+".pe.bam").getAbsolutePath()).append(" &");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            
            BufferedReader br = IOUtils.getTextReader(bwaScriptS);
            System.out.println(br.readLine());
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        /**************************** sort **************************************/
        try{
            BufferedWriter bw = IOUtils.getTextWriter(sortScriptS);
            for (int i = 0; i < names.length; i++) {
       //samtools sort -m 1G -o /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.pos.bam -O bam -@ 10 /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.bam &
                StringBuilder sb = new StringBuilder("samtools sort -m ");
                sb.append(sortMemory).append(" -o ").append(bamfileDirS).append(names[i]).append(".sorted.bam").append(" -O bam -@ ").append(sortThread).append(" ").append(bamfileDirS).append(names[i]).append(".pe.bam && ");
                sb.append("samtools index ").append(bamfileDirS).append(names[i]).append(".sorted.bam &");
                bw.write(sb.toString());bw.newLine();
            }
            bw.flush();bw.close();
            
            BufferedReader br = IOUtils.getTextReader(sortScriptS);
            System.out.println(br.readLine());
            br.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        /**************************** merge **************************************/
        try{
            BufferedWriter bw = IOUtils.getTextWriter(mergeScriptS);
            
            //samtools merge out.bam in1.bam in2.bam in3.bam
            StringBuilder sb = new StringBuilder("samtools merge ");
            sb.append(bamfileDirS).append(mergeName).append(" ");
            for (int i = 0; i < names.length; i++) {
                sb.append(bamfileDirS).append(names[i]).append(".sorted.bam").append(" ");
            }
            sb.deleteCharAt(sb.length()-1);
            sb.append(" && mv ").append(bamfileDirS).append(mergeName).append(" ").append(mergefileDirS);
            sb.append(" && ");
            sb.append("samtools index ").append(mergefileDirS).append(mergeName);
            bw.write(sb.toString());
            bw.newLine();
            bw.flush();bw.close();
            
            BufferedReader br = IOUtils.getTextReader(mergeScriptS);
            System.out.println(br.readLine());
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        /**************************** bam QC **************************************/
        try{
            BufferedWriter bw = IOUtils.getTextWriter(bamQCScriptS);
            //samtools stats -r /data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa 180803.pe.sorted.bam > 180803.pe.sorted.bam.bc
            StringBuilder sb = new StringBuilder("samtools stats -r ");
            sb.append(RefSeqS).append(" ").append(mergefileDirS).append(mergeName).append(" > ");
            sb.append(mergefileDirS).append(mergeName).append(".bc");
            sb.append(" && ").append("mv ").append(mergefileDirS).append(mergeName).append(".bc ").append(statfileDirS);
            sb.append(" && ");
            sb.append("perl /data1/programs/samtools-1.8/misc/plot-bamstats -p ").append(statfileDirS).append(mergeName).append(".bc ")
                    .append(statfileDirS).append(mergeName).append(".bc");
            bw.write(sb.toString());
            bw.newLine();
            
            sb = new StringBuilder();
            sb.append("samtools flagstat ").append(mergefileDirS).append(mergeName).append(" > ").append(mergefileDirS).append(mergeName).append(".flagstat.txt &");
            bw.write(sb.toString());
            bw.newLine();
            
            
            bw.flush();bw.close();
            
            BufferedReader br = IOUtils.getTextReader(bamQCScriptS);
            System.out.println(br.readLine());
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    
    public void samtoolsMerge(){
        //String inputDirS = "/Users/Aoyue/Downloads/huadagene/bamfile";
        //String outputPerlS = "/Users/Aoyue/Downloads/huadagene/runMerge.sh";
        String inputDirS = "/data1/home/aoyue/huada/bamfile";
        String outputPerlS = "/data1/home/aoyue/huada/runMerge.sh";
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".bam");
        Arrays.sort(fs);
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outputPerlS);
            bw.write("samtools merge mergeWheat24SM.bam ");
            StringBuilder sb = new StringBuilder();
            for(int i = 0; i < fs.length; i++){
                String inBam = fs[i].getAbsolutePath();
                System.out.println(inBam);
                sb.append(inBam).append(" ");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.write(" &");
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void samtoolsSort(){
        //samtools sort -m 6G -o 180803.pe.sorted.bam -O bam -@ 8 180803.pe.bam && echo "** BAM sort done **"
        String inputDirS = "/data1/home/aoyue/huada/bamfile/";
        String outputPerlS = "/data1/home/aoyue/huada/runSort.pl";
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".bam");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            String fullName = fs[i].getName();
            String[] all = fullName.split(".pe");
            nameSet.add(all[0]);
        }
        String[] names = nameSet.toArray(new String[nameSet.size()]);
        Arrays.sort(names);
        int memory = 8;
        int numthreads = 8;
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outputPerlS);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder("samtools sort");
                sb.append(" -m ").append(memory).append("G").append(" -o ");
                sb.append(new File(inputDirS,names[i] + ".pe.sorted.bam" ).getAbsolutePath());
                sb.append(" -O ").append("bam").append(" -@ ").append(numthreads).append(" ");
                sb.append(new File (inputDirS,names[i] + ".pe.bam" ).getAbsolutePath());
                sb.append(" && ");
                sb.append("echo \"** BAM sort done **\"");
                sb.append(" && ");
                sb.append("samtools index ");
                sb.append(new File(inputDirS,names[i] + ".pe.sorted.bam" ).getAbsolutePath());
                sb.append(" && ");
                sb.append("echo \"** BAM index done **\"").append(" &");
                String cmd = sb.toString();
                //String command = "system(\""+cmd+"\");";
                String command = cmd;
                bw.write(command);
                bw.newLine();
            }
            bw.flush();
            bw.close();
//            StringBuilder sb = new StringBuilder("perl ");
//            sb.append(outputPerlS);
//            Runtime run = Runtime.getRuntime();
//            Process p = run.exec(sb.toString());
//            System.out.println(sb.toString());
//            BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//            p.waitFor();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void testsplit(){
        String inDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "1.clean.fq.gz");
        String outfileS = "/Users/Aoyue/Documents/list.sh";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < fs.length; i++) {
                sb = new StringBuilder();
                String fullName = fs[i].getName();
                String[] all = fullName.split("_",4);
                String sample = all[3];
                sb.append(fs[i].getName());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(fs.length);  
    }
    
    private void listSpecificalFiles() {
        String inDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "fq.gz");
        String outfileS = "/Users/Aoyue/Documents/Datalist.txt";
        //String[] header = {"FileName", "FileSize"};
        String[] header = {"FileName"};
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < header.length; i++) {
                sb.append(header[i]).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                sb = new StringBuilder();
                
                double qq = (double)fs[i].length()/1024/1024/1024;
                DecimalFormat df = new DecimalFormat ("0.00");

                sb.append(fs[i].getName());
                        //.append("\t").append(df.format(qq));
//                sb.append(subFs[i].getName()).append("\t").append(subFs[i].getAbsolutePath()).append("\t").append((double)subFs[i].length()/1024/1024/1024).append("\t").append("350bp");
                
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(fs.length);    
    }
    
    private void alignBWA () {
        String infileS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/001_script/SM8.t.txt";
        String thread = "16";  //根据实际情况修改，使用的线程数
        String indexFileS = "/data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz";
        String inputDirS = "/data2/aoyue/a-wheatRawdata1219";
        String outputDirS = "/data2/aoyue/a-output/bamfile/";
        String ScriptS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/001_script/001_bwa.sh";
        
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        Collections.sort(namelist);
        String[] names = namelist.toArray(new String[namelist.size()]);
        
        //bwa mem -t 100 -R '@RG\tID:HRV-L1\tPL:illumina\tSM:HRV-L1\tLB:HRV-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD/HRV-L1_1.fq.gz /data2/sharedData/Jiao/ABD/HRV-L1_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HM.pe.bam && echo "** bwa mapping done **" &
        try {
            BufferedWriter bw = IOUtils.getTextWriter(ScriptS);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder("bwa mem -t ");
                sb.append(thread).append(" -R ").append("'@RG\\tID:").append(names[i]).append("\\t").append("PL:illumina").append("\\t").append("SM:").append(names[i]).append("\\t").append("LB:").append(names[i]).append("' ");
                sb.append(indexFileS).append(" ");
                sb.append(new File(inputDirS, names[i]+"_1.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(inputDirS, names[i]+"_2.fq.gz").getAbsolutePath());
                sb.append(" | ").append("samtools view -S -b - > ").append(new File(outputDirS, names[i]+".pe.bam").getAbsolutePath()).append(" &");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        try{
            BufferedReader br = IOUtils.getTextReader(ScriptS);
            System.out.println(br.readLine());
            br.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void fastQC () {
//        String inputDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/";
//        String outputDirS = "/Users/Aoyue/Downloads/huadagene/fastqc/";
        String inputDirS = "/Users/Aoyue/Documents/P18Z12200N0030_WHEsikR/Clean/WHYD18111796_A/";
        String outputDirS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/000_fastqc";
        
        try {
            StringBuilder sb = new StringBuilder("/Users/Aoyue/software/FastQC/fastqc");
            File[] fs = new File (inputDirS).listFiles();
            fs =  IOUtils.listFilesEndsWith(fs, ".gz");
            for (int i = 0; i < fs.length; i++) {
                sb.append(" ").append(fs[i].getAbsoluteFile());
            }
            sb.append(" -o ").append(outputDirS);
            String cmd = sb.toString();
            System.out.println(cmd);
            Runtime run = Runtime.getRuntime();
            Process p = run.exec(cmd);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
            System.out.println("Fastqc evalutation is finished at" + outputDirS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    private void checkMd5() {
        //ori 文件 54de998c7883a4592a1683bec2590d64  K16HL0119_1_clean.fq.gz
        //ori 文件 0dc05f49a968937197f20841ef8427b2  Clean/WHYD18111796_A/V300010959_L1_B5RDWHEsikRAAAAA-533_1.fq.gz
        //des文件 MD5 (/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/180803_I13_V100004234_L1_WHEkapRAADT-585_1.clean.fq.gz) = 472e1633b40b0ae3866a820f5b8637a5
        //des文件 MD5 (/Users/Aoyue/Documents/P18Z12200N0030_WHEsikR/Clean/WHYD18111796_A/V300010959_L1_B5RDWHEsikRAAAAA-533_1.fq.gz) = 0dc05f49a968937197f20841ef8427b2
        String ori = "/Users/Aoyue/Documents/P18Z12200N0030_WHEsikR/md5.txt"; //原始md5
        String des = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/Wheat20181219.md5.txt"; //mac生成的md5文件
        HashMap<String, String> fmd5Map = new HashMap<>(); //建立一个键值对应的hashmap,此时hashmap为空。下文会把原始的ori文件放入hashmap中去
        TableInterface oT = new RowTable(ori, "  "); //分隔符是2个空格，将ori文件读进表格
        TableInterface dT = new RowTable(des, " "); //分隔符是1个空格，将des文件读进表格
        String choice = "YES";
        try{
            BufferedReader br = IOUtils.getTextReader(ori); //先读进原始生成的md5文件
            String temp;
            List<String> l = null;
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp, "  ");
                //fmd5Map.put(l.get(1), l.get(0));
                fmd5Map.put(l.get(1).split("_A/")[1], l.get(0));
            }
            br.close();
            int cnt = 0;
            br = IOUtils.getTextReader(des); //读进mac生成md5文件
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp, " ");
                //String key = l.get(1).split("00010496/")[1].replaceFirst("\\)", "");
                String key = l.get(1).split("_A/")[1].replaceFirst("\\)", "");
                String value = fmd5Map.get(key);
                /********* 文件全部下载完毕，可以进行一一对应 ***********/
                if (value == null) { //如果value为空，则为真。
                System.out.println(key+"\tdoesn't exist");
                continue;
                }
                if (value.equals(l.get(3))) {
                    cnt++;
                    continue;
                }
                System.out.println(key + "\t is incorrect");
            }
            System.out.println(cnt + " key is correct");
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void mkMd5(){
        /*该方法只用将文件输入路径和输出文件名进行修改，即可使用*/
//        String inputDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496";
//        String outfileS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/checkmd5.txt";
        //String inputDirS = "/Users/Aoyue/Documents/Data/project/WheatGBSVI/02IlluminaSeq/data/2.cleandata/20180601-51-NEB12_TKD180600155";
        //String inputDirS = "/Users/Aoyue/Documents/P18Z12200N0030_WHEsikR/Clean/WHYD18111796_A";
        String inputDirS = "/Volumes/LuLab3T_42/a-wheatRawdata1219";
        String outfileS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/Wheat20181219.md5.txt";
        String suffix = ".gz";
        
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            /*1#新建存放数据的文件对象fall，将该目录下所有以fq.gz结尾的文件名列出来，并存放在fs文件数组中*/
            File fall = new File (inputDirS);
            File[] fs = IOUtils.listRecursiveFiles(fall);
            fs =  IOUtils.listFilesEndsWith(fs, suffix);
            /*2# 根据fq.gz结尾的文件，生成md5计算的脚本，输入到StringBuilder sb中（这里用循环的方法）*/
            for (int i = 0; i < fs.length; i++) {
                sb.append("md5").append(" ").append(fs[i].getAbsoluteFile()).append("\n");
            }
            /*3# 将生成的sb脚本转换成String类型，对象名为cmd，并打印出来检查是否无误*/
            String cmd = sb.toString();
            System.out.println(cmd);
            /*4# 调用终端命令，运行md5命令，并将输出的结果写入输出文件outfileS中*/
            Runtime run = Runtime.getRuntime();
            Process p = run.exec(cmd);           
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream())); 
            //p.getInputStream()返回InputStream类型；InputStreamReader构造器函数需要加 InputStream 类型的参数； 
            //InputStreamReader类继承了Reader类. BufferedReader构造器需要加Reader类的参数。
            StringBuilder ssb = new StringBuilder();
            String line;
            while ((line = br.readLine()) != null) {
            ssb.append(line).append("\n");
            }
            String result = ssb.toString();
            System.out.println(result);  
            p.waitFor();
            bw.write(result);
            bw.flush();
            bw.close();            
            System.out.println("md5Check is finished at " + outfileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);       
        }       
    }
    
    
    public static void main (String[] args){
        new WheatBWA();
        System.out.println("\n**********************************" );
        System.out.println("Here is the main class of WheatBWA" );
        
    }  
}

