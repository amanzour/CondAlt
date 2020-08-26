/*
 * This program reads a *.plot file and constructs a vienna structure
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package unafold2vienna;

import java.io.BufferedReader;
import java.io.FileReader;

/**
 *
 * @author manzouro
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String seq="";
        try{
        BufferedReader in = new BufferedReader(new FileReader(args[0]));
        String str;
        while ((str = in.readLine()) != null)
            if (!str.startsWith(">"))
                seq=str;
        }
        catch (Exception e){
              System.err.println("Error: Fasta File Problem" + e.getMessage());
           }
        int L=seq.length();
        char[] repeat = new char[L];
         for(int n=0; n<L; n++)
            repeat[n]='.'; 
        try{
        BufferedReader in = new BufferedReader(new FileReader(args[1]));
        String str;
        int i=0;int j=0;
        while ((str = in.readLine()) != null)
            if (!str.startsWith("i")){
               i=Integer.parseInt(str.split("\t")[0]);
               j=Integer.parseInt(str.split("\t")[1]);
               repeat[i-1]='(';repeat[j-1]=')';
            }
        }
        catch (Exception e){
              System.err.println("Error: plot File Problem" + e.getMessage());
           }
        String structure = new String(repeat);
        System.out.println(structure);
    }
}