/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGBS;

import java.io.IOException;

/**
 *
 * @author Aoyue
 */
public class GBSEntrance {
    public GBSEntrance() throws IOException{
        //this.PlateAndID();
        this.DataProcessor();
        
        
    }
    public void PlateAndID(){
        new PlateAndID();
        
    }
    public void DataProcessor() throws IOException{
        new DataProcessor();
        //new GBScp();
    }
    public static void main (String[] args) throws IOException{
        new GBSEntrance();
        
    }  
    
}
