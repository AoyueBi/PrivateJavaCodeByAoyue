/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGBS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class DataProcessor {
    public DataProcessor() throws IOException{
        //this.sample();
       //this.checkMd51();
        //this.runtime();
       this.checkMd5();
       this.mergeFile();
        
    }
    
    public void checkMd5(){
        /*该方法只用将文件输入路径和输出文件名进行修改，即可使用*/
        String inputDirS = "/Users/Aoyue/Documents/testCheckMd5/";
        String outfileS = "/Users/Aoyue/Documents/checkmd52.txt";
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
            System.out.println("md5Check is finished at" + outfileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);       
        }       
    }
    public void mergeFile(){
        String inputDirS = "/Users/Aoyue/Documents/testCheckMd5/";
        String outfileS = "/Users/Aoyue/Documents/orimd5_merge.txt";
        File fall = new File (inputDirS);
        File[] orimd5 = IOUtils.listRecursiveFiles(fall);
        orimd5 = IOUtils.listFilesEndsWith(orimd5, ".md5");
        List<File> orimd5List = Arrays.asList(orimd5);
        try{
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for(int i = 0; i< orimd5.length; i++){
                BufferedReader br = IOUtils.getTextReader(orimd5List.get(i).getPath());
                String temp = null;
                while((temp = br.readLine())!= null){
                    sb.append(temp).append("\n");                   
                }                
                String a = orimd5List.get(i).getName();
                System.out.println(a);
            }
            bw.write(sb.toString());
            bw.flush();
            bw.close();
            //br.close();            
        }
        catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }       
    }
    
    public void runtime() throws IOException{         
        try {
            StringBuilder sb = new StringBuilder("/sbin/md5");
            BufferedWriter bw = IOUtils.getTextWriter("/Users/Aoyue/Documents/checkMd5.txt");
            sb.append(" ").append("/Users/Aoyue/Documents/testCheckMd5/data1_00/SRR1575344_1.fq.gz").append(" ").append(">").append(" ").append("/Users/Aoyue/Documents/testCheckMd5.txt");
            String cmd = sb.toString();
            System.out.println(cmd);
            Runtime run = Runtime.getRuntime();   
            Process p = run.exec(cmd);           
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            StringBuffer ssb = new StringBuffer();
            String line;
            while ((line = br.readLine()) != null) {
            ssb.append(line).append("\n");
            }
            String result = ssb.toString();
            bw.write(result);
            bw.flush();
            bw.close();
          //System.out.println(result); 
            p.waitFor();
        }
        catch (Exception e) {
        System.out.println(e.toString());
        System.exit(1);
        }       
        System.out.println("md5Check is finished");
    }
    
    public void checkMd51(){
        String infileDirS = "/Users/Aoyue/Documents/testCheckMd5";
        //String outfileS = "/Users/Aoyue/Documents/testCheckMd5.txt";
        File fsall = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(fsall);
        fs = IOUtils.listFilesEndsWith(fs, "fq.gz");
        List<File> fsList = Arrays.asList(fs);
        for(int i = 0; i< fs.length; i++){
            String a = fsList.get(i).getName();
            System.out.println(a);
        }  
        fsList.parallelStream().forEach(f -> {
            StringBuilder sb = new StringBuilder("/sbin/md5");
            sb.append(" ").append(f.getPath());
            String cmd = sb.toString();
            System.out.println(cmd);       
            Runtime run = Runtime.getRuntime();        
            try {
                Process p = run.exec(cmd);
                p.waitFor();
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                StringBuffer ssb = new StringBuffer();
                String line;
                while ((line = br.readLine()) != null) {
                ssb.append(line);
                }
                String result = ssb.toString();
                System.out.println(result);              
            } 
            catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }            
        });        
        System.out.println("md5Check is finished");       
    }
   
    public void sample(){
        String infileS = "/Users/Aoyue/Downloads/2.cleandata/20180601-51-NEB12_TKD180600155/20180601-51-NEB12_TKD180600155_1.clean.fq.gz";
        String outfileS = "/Users/Aoyue/Downloads/20180601-51-NEB12_TKD180600155_1.clean.fq";
        int length = 120000;
        try{
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < length; i++){
                bw.write(br.readLine());
                bw.newLine();
            } 
            bw.flush();
            bw.close();
            br.close();          
        }
        catch(Exception e){
            e.printStackTrace();        
        }  
    }
}