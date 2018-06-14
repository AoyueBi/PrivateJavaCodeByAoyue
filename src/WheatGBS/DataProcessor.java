/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGBS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class DataProcessor {
    public DataProcessor(){
        this.sample();
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