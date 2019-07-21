/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import aoyue.utils.Monitor;
import java.io.File;
import vcf.VCFtools;

/**
 *
 * @author Aoyue
 */
public class MainEntrance {
    
    public MainEntrance(){
        new VCFtools();
        
        
    }
    
    public static void main (String[] args) throws InterruptedException{
        //System.out.println("******************************************************" );
        //System.out.println("Here is the main class of Basic Genetics Analysis !\n" );
        
        
        System.out.println("\nHello! To use this monitor jar, please add 4 parameters and prepare two files.\n");
        System.out.println("@param1: keepjobnum eg:20 \n@param2: sleeptime eg:5 \n@param3: getJobRunshS eg:/Users/Aoyue/Documents/getJobRun.sh \n@param4: allshS eg:/Users/Aoyue/Documents/001_allsh_samtoolsView.sh\n");
        System.out.println("From linux command line: java -jar Aoyue-monitor.jar 20 5 /Users/Aoyue/Documents/getJobRun.sh /Users/Aoyue/Documents/001_allsh_samtoolsView.sh > log_monitor.txt & \n");
        System.out.println("Note:\n@param1: The job you wanna keep in linux process \n@param2: Unit: min, the interval time you"
        + " wanna check the ps aux process\n@param3: The sh file to produce the job num running now eg:ps aux|grep aoyue|grep samtools|grep view|wc -l \n"
        + "@param4: All the CMDs in one script, one line one sample, please DO NOT add & symbol at the end of each line.\n");

        new Monitor(args[0],args[1],args[2],args[3]);
        //new Monitor().monitor("6", "6", "fbgb", "shuch");
                
       //new MainEntrance();
        
    }
    
}
