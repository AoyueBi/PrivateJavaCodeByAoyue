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
//        this.checkMd5();
//        this.fastQC();
        this.alignBWA();
        //this.listSpecificalFiles();
        //this.testspilt();
        //this.samtoolsSort();
        //this.samtoolsMerge();
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
        String samtoolsPath = "samtools";
        //String inputDirS = "/Users/Aoyue/Downloads/huadagene/bamfile";
        //String outputPerlS = "/Users/Aoyue/Downloads/huadagene/runsorted.pl";
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
        //bwa mem -t  nthreads ref.fa read1.fq read2.fq > aln-pe.sam
        /*bwa mem  -t 8 -R '@RG\tID:foo\tPL:illumina\tSM:WHB5EXONPEP00010496\tLB:LB' 
        /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data1/home/aoyue/testSampleData/180803_1.fq.gz 
        /data1/home/aoyue/testSampleData/180803_2.fq.gz | samtools view -S -b - > /data1/home/aoyue/huada/bamfile/180803.pe.bam 
        && echo "** bwa mapping done **" &
        */
        //int numCores = Runtime.getRuntime().availableProcessors();
        int memory = 8;
        int numthreads = 8;
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outputPerlS);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder(samtoolsPath+" sort");
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
    
    public void testspilt(){
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
        String bwaPath = "bwa";
        //String indexFileS = "/Users/Aoyue/Documents/Data/wheat/ABD/bwaLib/abd_iwgscV1.fa.gz";
        String indexFileS = "/data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz";
        //String inputDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496";
        String inputDirS = "/data1/home/aoyue/WHB5EXONPEP00010496";
        //String outputDirS = "/Users/Aoyue/Downloads/huadagene/sam/";
        String outputDirS = "/data1/home/aoyue/huada/bamfile/";
        String outputPerlS = "/Users/Aoyue/Downloads/huadagene/runAlign.pl";
        //String outputPerlS = "/data1/home/aoyue/huada/runAlign.pl";
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "1.clean.fq.gz");
        
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            String fullName = fs[i].getName();
            String[] all = fullName.split("_1");
            nameSet.add(all[0]);
        }
        String[] names = nameSet.toArray(new String[nameSet.size()]);
        Arrays.sort(names);
        for (int i=0; i< fs.length;i++){
            System.out.println(names[i]);
        }
        //bwa mem -t  nthreads ref.fa read1.fq read2.fq > aln-pe.sam
        /*bwa mem  -t 8 -R '@RG\tID:foo\tPL:illumina\tSM:WHB5EXONPEP00010496\tLB:LB' 
        /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data1/home/aoyue/testSampleData/180803_1.fq.gz 
        /data1/home/aoyue/testSampleData/180803_2.fq.gz | samtools view -S -b - > /data1/home/aoyue/huada/bamfile/180803.pe.bam 
        && echo "** bwa mapping done **" &
        */
        int numCores = Runtime.getRuntime().availableProcessors();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outputPerlS);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder(bwaPath+" mem");
                sb.append(" -t ").append(numCores).append(" -R ").append("'@RG\\tID:").append(names[i]).append("\\t")
                        .append("PL:illumina").append("\\t").append("SM:").append(names[i]).append("\\t").append("LB:").append(names[i]).append("' ")
                        .append(indexFileS).append(" ");
                sb.append(new File(inputDirS, names[i]+"_1.clean.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(inputDirS, names[i]+"_2.clean.fq.gz").getAbsolutePath());
                sb.append(" | ").append("samtools view -S -b - > ").append(new File(outputDirS, names[i]+".pe.bam").getAbsolutePath()).append(" &");
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
    
    private void fastQC () {
        String inputDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/";
        String outputDirS = "/Users/Aoyue/Downloads/huadagene/fastqc/";
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
        //des文件 MD5 (/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/180803_I13_V100004234_L1_WHEkapRAADT-585_1.clean.fq.gz) = 472e1633b40b0ae3866a820f5b8637a5
        String ori = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/md5.txt"; //原始md5
        String des = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/checkmd5.txt"; //mac生成的md5文件
        HashMap<String, String> fmd5Map = new HashMap<>(); //建立一个键值对应的hashmap,此时hashmap为空。下文会把原始的ori文件放入hashmap中去
        TableInterface oT = new RowTable(ori, "  "); //分隔符是2个空格，将ori文件读进表格
        TableInterface dT = new RowTable(des, " "); //分隔符是1个空格，将des文件读进表格
        try{
            BufferedReader br = IOUtils.getTextReader(ori);
            String temp;
            List<String> l = null;
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp, "  ");
                fmd5Map.put(l.get(1), l.get(0));
            }
            br.close();
            int cnt = 0;
            br = IOUtils.getTextReader(des);
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp, " ");
                String key = l.get(1).split("00010496/")[1].replaceFirst("\\)", "");
                String value = fmd5Map.get(key);
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
        String inputDirS = "/Volumes/AoyueKe/maizeAGPv4_GERP";
        String outfileS = "/Volumes/AoyueKe/maizeAGPv4_GERP/md5.txt";
        
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            /*1#新建存放数据的文件对象fall，将该目录下所有以fq.gz结尾的文件名列出来，并存放在fs文件数组中*/
            File fall = new File (inputDirS);
            File[] fs = IOUtils.listRecursiveFiles(fall);
            fs =  IOUtils.listFilesEndsWith(fs, ".gz");
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
        
    }  
}

