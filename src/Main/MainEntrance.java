/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import vcf.VCFtools;

/**
 *
 * @author Aoyue
 */
public class MainEntrance {
    
    public MainEntrance(){
        new VCFtools();
        
        
    }
    
    public static void main (String[] args){
        System.out.println("******************************************************" );
        System.out.println("Here is the main class of Basic Genetics Analysis !" );
        new MainEntrance();
        
    }
    
}
