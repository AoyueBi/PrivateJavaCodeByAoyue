/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Monitor {
    
    String keepjobnum =null;
    String sleeptime =null;
    String getJobRunshS=null;
    String allshS =null;
    
    public Monitor(String jobnum, String sleepmin, String getJobRunshS,String allshS) throws InterruptedException{
        this.monitor(jobnum,sleepmin,getJobRunshS,allshS);
        
    }
    public Monitor(){
        
    }
    
    public void monitor(String keepjobnum, String sleeptime, String getJobRunshS,String allshS) throws InterruptedException{
        /************* parameters and file ********/
        //********parameters
        //int jobnum = 4; //根据脚本里的线程数，设置cluster每次提交的任务数,保证任务数一直是 jobnum
        //int sleepmin = 5; //根据一个样本运行的完成时间来定
        //String scriptDirS = null; //
        //********file
        //String allshS = null;  //所有待运行的命令写进一个脚本里，注意每行不要使用 "&"符号；
        //String getJobRunshS = null; //是一个脚本，包含一行命令，用来获取cluster正在执行的任务数，从而得到一个int类型的数字。eg: ps aux|grep aoyue|grep samtools|grep view|wc -l
        
        /********159.226.116.204 test********/
//        allshS = "/data4/home/aoyue/vmap2/ab/test_monitor/001_bamdata/001_samtoolsView.sh";
//        getJobRunshS = "/data4/home/aoyue/vmap2/ab/test_monitor/001_bamdata/getJobRun.sh";
//        scriptDirS = "/data4/home/aoyue/vmap2/ab/test_monitor/001_bamdata/scriptDirS";

        int jobnum = Integer.parseInt(keepjobnum);
        int sleepmin = Integer.parseInt(sleeptime);
        String scriptDirS = new File(new File (allshS).getParent(),"scriptDirS").getAbsolutePath();
        new File(scriptDirS).mkdirs();
        
        try{
            /*#########################获取所有命令，存入HashMap中##############################*/
            /*#########################获取文件行数，得到总tjobs################################*/
            /**
             * pseudo-code:建立一个HashMap,第一行key是0001，value是第一行的代码；依次类推；
             */
            BufferedReader brall = IOUtils.getTextReader(allshS);
            HashMap<String,String> hm = new HashMap<>();
            String temp;
            int cnt=0;
            while((temp =brall.readLine()) != null){
                cnt++;
                String k = PStringUtils.getNDigitNumber(4, cnt);
                hm.put(k, temp);
            }
            brall.close();
            int tjobs = cnt;
            System.out.println("***********************************************");
            System.out.println(new SimpleDateFormat().format(new Date()));
            System.out.println("There are  " + tjobs + "  total jobs");
            System.out.println("================================================" + "\n");
            int remainJob =tjobs;

            /*#########################获取linux系统正在运行的任务数，即linux job数##############################*/
            
            for(;;){
                TimeUnit.MINUTES.sleep(1); //让程序休息1分钟，等待命令提交后再进行ps监控；
                String cmd = "sh " + getJobRunshS;
                Process p = Runtime.getRuntime().exec(cmd);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                int runningJob = Integer.parseInt(br.readLine());
                br.close();
                p.waitFor();
                
                
                //System.out.println("***********************************************");
                //System.out.println(new SimpleDateFormat().format(new Date()));
                //System.out.println("There are " + runningJob + "   jobs running in ps process");
                //System.out.println("================================================" + "\n");
                /**
                 * 伪代码：集群保持一直有5个命令在运行，当有4个命令在运行时，添加1个到linux； 当有2个命令在运行时，自动添加3个命令到linux；
                 * 当有5个命令在运行时，不添加任何命令，程序休眠，并于5分钟后再次进行check。
                 * if(runningJob < 5)
                 */
                
                /*######################### 开始提交任务并运行 ###################*/
                if(runningJob < jobnum){
                    int submitjob = jobnum - runningJob; 
                    //接下来对应该提交任务数 和 剩余任务数进行比较.最开始，剩余任务书等于脚本内的命令数，后来每读一行，就少一个任务。
                    //if(提交任务数小于剩余任务数，就提交 提交的任务数； 否则提交 剩余任务数)
                    if(submitjob <= remainJob){
                        /******* 多线程提交任务 *********/
                        BufferedWriter[] bw = new BufferedWriter[submitjob];
                        File[] fs = new File[submitjob];
                        for(int i=0;i<submitjob;i++){
                            String linenumber = PStringUtils.getNDigitNumber(4, 1+tjobs-remainJob);
                            String scriptname = linenumber + ".sh";
                            fs[i]=new File(scriptDirS,scriptname);
                            bw[i] = IOUtils.getTextWriter(fs[i].getAbsolutePath());
                            String indicmd = hm.get(linenumber);
                            bw[i].write(indicmd + "\n");bw[i].flush();bw[i].close();
                            remainJob--; 
                        }
                        List<File> fsList = Arrays.asList(fs);
                        fsList.parallelStream().forEach(f -> {
                            String CMD = "sh " + f.getAbsolutePath() + " &";
                            System.out.println(new SimpleDateFormat().format(new Date()) + "    The " + f.getName() + " has been submitted \n"+ CMD);
                            Process ps = null;
                            try {
                                ps = Runtime.getRuntime().exec(CMD);
                            } catch (IOException ex) {
                                Logger.getLogger(Monitor.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        });
                    }
                    else{//如果提交的命令数 大于剩余的命令数，比如还有3个要提交，但是脚本里只剩下1个命令没运行了，此时我们只提交最后剩余的几个；
                        BufferedWriter[] bw = new BufferedWriter[remainJob];
                        File[] fs = new File[remainJob];
                        int circle = remainJob;
                        for(int i=0;i<circle;i++){ 
                            String linenumber = PStringUtils.getNDigitNumber(4, 1+tjobs-remainJob);
                            String scriptname = linenumber + ".sh";
                            fs[i]=new File(scriptDirS,scriptname);
                            bw[i] = IOUtils.getTextWriter(fs[i].getAbsolutePath());
                            String indicmd = hm.get(linenumber);
                            bw[i].write(indicmd + "\n");bw[i].flush();bw[i].close();
                            remainJob--; 
                        }
                        List<File> fsList = Arrays.asList(fs);
                        fsList.parallelStream().forEach(f -> {
                            String CMD = "sh " + f.getAbsolutePath() + " &";
                            System.out.println(new SimpleDateFormat().format(new Date()) + "    The " + f.getName() + " has been submitted \n"+ CMD);
                            Process pc = null;
                            try {
                                pc = Runtime.getRuntime().exec(CMD);
                            } catch (IOException ex) {
                                Logger.getLogger(Monitor.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        });
                    }
                }
                else if(runningJob >= jobnum){ //提交任务数大于jobnum规定数目时
                    //System.out.println("SSSSSSSSS");
                    //System.out.println("Now I will be sleeping for"+ sleepmin + "min\n");
                    TimeUnit.MINUTES.sleep(sleepmin); //让程序休眠5分钟
                }
                if(remainJob ==0) {
                    System.out.println("***********************************************");
                    System.out.println(new SimpleDateFormat().format(new Date()));
                    System.out.println("The remainJob is 0, this monitor process will quit. Good luck!");
                    System.out.println("================================================" + "\n");
                    break;
                }
            }//无限循环
        }
        catch(IOException ex){
            ex.printStackTrace();
            System.exit(1);
        }
    }
}
