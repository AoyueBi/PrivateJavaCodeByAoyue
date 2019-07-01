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
    public Monitor(){
        
    }
    public void monitor() throws InterruptedException{
        int jobnum = 4;
        /********local test********/
//        String allshS="/Users/Aoyue/Documents/001_samtoolsView.sh";
//        String getJobRun = "/Users/Aoyue/Documents/testlinux.sh";
//        String jobRun = "/Users/Aoyue/Documents/jobRun.txt";
        /********159.226.116.204 test********/
        String allshS = "/data4/home/aoyue/vmap2/ab/test_monitor/001_bamdata/001_samtoolsView.sh";
        String getJobRun = "/data4/home/aoyue/vmap2/ab/test_monitor/001_bamdata/getJobRun.sh";
        //String jobRun = "/data4/home/aoyue/vmap2/ab/test_monitor/001_bamdata/jobRun.txt";
        String scriptDirS = "/data4/home/aoyue/vmap2/ab/test_monitor/001_bamdata/scriptDirS";
        new File(scriptDirS).mkdirs();
        
        try{
            /*#########################获取文件行数，即总共需要运行的job数##############################*/
            LineNumberReader reader = new LineNumberReader(new FileReader(new File(allshS)));
            reader.skip(Long.MAX_VALUE);
            int tjobs = reader.getLineNumber();
            reader.close();
            System.out.println(tjobs + "    total jobs");
            System.out.println("================================================" + "\n");
            int remainJob =tjobs;
            
            /*#########################获取所有命令，存入HashMap中##############################*/
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
            /*#########################获取linux系统正在运行的任务数，即linux job数##############################*/
            
            for(;;){
                //BufferedWriter bw = IOUtils.getTextWriter(jobRun);
                String cmd = "sh " + getJobRun;
                Process p = Runtime.getRuntime().exec(cmd);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                int runningJob = Integer.parseInt(br.readLine());
                //bw.write(runningJob + "\n");bw.flush();bw.close();
                br.close();
                p.waitFor();
                System.out.println("***********************************************");
                System.out.println(new SimpleDateFormat().format(new Date()));
                System.out.println(runningJob + "   jobs are running in ps process");
                System.out.println("================================================");
                System.out.println("================================================" + "\n");
                /**
                 * 伪代码：集群保持一直有5个命令在运行，当有4个命令在运行时，添加1个到linux； 当有2个命令在运行时，添加3个命令到linux；
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
                            System.out.println("The " + f.getName() + " is \n"+ CMD);
                            //System.out.println(new SimpleDateFormat().format(new Date()));
                            //System.out.println("================================================" + "\n");
                            Process ps = null;
                            try {
                                ps = Runtime.getRuntime().exec(CMD);
                            } catch (IOException ex) {
                                Logger.getLogger(Monitor.class.getName()).log(Level.SEVERE, null, ex);
                            }
//                            try {
//                                ps.waitFor();
//                            } catch (InterruptedException ex) {
//                                Logger.getLogger(Monitor.class.getName()).log(Level.SEVERE, null, ex);
//                            }
                            
                        });
                        
                        /******* 此方法为单线程，故需要设置为多线程运行 *********/
//                        for(int i=0;i<submitjob;i++){ //开始读入文件，将文件写入脚本中并执行！
//                            String scriptname = PStringUtils.getNDigitNumber(4, 1+tjobs-remainJob) + ".sh";
//                            bw[i] = IOUtils.getTextWriter(new File(scriptDirS,scriptname).getAbsolutePath());
//                            bw[i].write(br.readLine() + "\n");bw[i].flush();bw[i].close();
//                            remainJob--; 
//                            cmd = "sh " + new File(scriptDirS,scriptname).getAbsolutePath() + " &"; //注意命令里不要用后台运行号，在脚本里用后台运行号
//                            System.out.println("The " + scriptname + " is \n"+ cmd);
//                            System.out.println(new SimpleDateFormat().format(new Date()));
//                            System.out.println("The remaining jobs are " + remainJob);
//                            System.out.println("================================================" + "\n");
//                            p = Runtime.getRuntime().exec(cmd);
//                            p.waitFor();
//                        }
                    }
                    else{//如果提交的命令数 大于剩余的命令数，比如还有3个要提交，但是脚本里只剩下1个命令没运行了，此时我们只提交最后剩余的几个；
                        /********多线程提交任务 ***********/
                        BufferedWriter[] bw = new BufferedWriter[remainJob];
                        File[] fs = new File[remainJob];
                        int circle = remainJob;
                        for(int i=0;i<circle;i++){ //
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
                            System.out.println("The " + f.getName() + " is \n"+ CMD);
                            //System.out.println(new SimpleDateFormat().format(new Date()));
                            //System.out.println("================================================" + "\n");
                            Process ps = null;
                            try {
                                ps = Runtime.getRuntime().exec(CMD);
                            } catch (IOException ex) {
                                Logger.getLogger(Monitor.class.getName()).log(Level.SEVERE, null, ex);
                            }
//                            try {
//                                ps.waitFor();
//                            } catch (InterruptedException ex) {
//                                Logger.getLogger(Monitor.class.getName()).log(Level.SEVERE, null, ex);
//                            }
                        });
                        
                            /******* 此方法为单线程，故需要设置为多线程运行 *********/
//                        BufferedWriter[] bw = new BufferedWriter[remainJob];
//                        for(int i=0;i<remainJob;i++){ //开始读入文件，将文件写入脚本中并执行！
//                            String scriptname = PStringUtils.getNDigitNumber(4, 1+tjobs-remainJob) + ".sh";
//                            bw[i] = IOUtils.getTextWriter(new File(scriptDirS,scriptname).getAbsolutePath());
//                            bw[i].write(br.readLine() + "\n");bw[i].flush();bw[i].close();
//                            remainJob--; 
//                            cmd = "sh " + new File(scriptDirS,scriptname).getAbsolutePath() + " &"; //注意命令里不要用后台运行号，在脚本里用后台运行号
//                            System.out.println("The " + scriptname + "is \n"+ cmd);
//                            System.out.println(new SimpleDateFormat().format(new Date()));
//                            System.out.println("The remaining jobs are " + remainJob);
//                            System.out.println("================================================" + "\n");
//                            p = Runtime.getRuntime().exec(cmd);
//                            p.waitFor();
//                        }
                    }
                }
                else if(runningJob >= jobnum){
                    //TimeUnit.MINUTES.sleep(5); //让程序休眠5分钟
                    TimeUnit.SECONDS.sleep(30);//让程序休眠30秒
                    System.out.println("**********************************");
                    System.out.println("Now I am sleeping for 30 seconds\n");
                    //TimeUnit.HOURS.sleep(1);//小时
                }
            }
            //无限循环
        }
        catch(IOException ex){
            ex.printStackTrace();
            System.exit(1);
        }
    }
}
