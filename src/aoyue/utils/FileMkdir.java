/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.utils;

import java.io.File;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class FileMkdir {
    public FileMkdir(){
        
    }
    public FileMkdir(String a, int b){
        //this.mkdir(a,b);
        
    }
    
    /**
     * 
     * @param ParentPath the path you want to mkdir
     * @param a the number you want to buid the file path
     */
    public void mkdir(String ParentPath, int a){
        //String outfileDirS = "/Users/Aoyue/Documents";
        for(int i =0 ; i < a; i++){
            String filename = PStringUtils.getNDigitNumber(3, i);
            File f = new File(ParentPath, filename);
            f.mkdirs();
        }
    }
}
