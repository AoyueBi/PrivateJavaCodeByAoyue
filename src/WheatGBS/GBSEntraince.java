/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGBS;

/**
 *
 * @author Aoyue
 */
public class GBSEntraince {
    public GBSEntraince(){
        //this.PlateAndID();
        this.DataProcessor();
        
        
    }
    public void PlateAndID(){
        new PlateAndID();
        
    }
    public void DataProcessor(){
        //new DataProcessor();
        new GBScp();
    }
    public static void main (String[] args){
        new GBSEntraince();
        
    }  
    
}
