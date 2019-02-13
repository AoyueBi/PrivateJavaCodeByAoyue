/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatUnclass;

/**
 *
 * @author Aoyue
 */
public class Test {
    int a = 4;
    
    public Test (int b) {
        this.a = b;
        System.out.println(this.a);
        System.out.println("#####");
        System.out.println(a);
    }
    
    public static void main (String[] args) {
        new Test(7);
    }
    
}