/*
 * This program reads a *.plot file and adds slippage = "slip" real number
 * which will be multiplied by length to become integer
 * It also has to read the length of the RNA as the third input
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package slip;

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
        int slipR=Integer.parseInt(args[1]);
        int L=Integer.parseInt(args[2]);
        //int slip = Math.round(slipR*L);
        int slip=0;
        if (slipR==-1)
            slip=(int)Math.ceil(1.0+Math.log10(L));
        else
            slip=slipR;
        try{
        BufferedReader in = new BufferedReader(new FileReader(args[0]));
        String str;
        int i=0;int j=0;
        while ((str = in.readLine()) != null)
            if ((!str.startsWith("i"))&&(Float.parseFloat(str.split("\t")[2])>0.5f)){
               i=Integer.parseInt(str.split("\t")[0]);
               j=Integer.parseInt(str.split("\t")[1]);
               System.out.println("P "+Integer.toString(i)+" "+Integer.toString(j)+" 1");
               for (int alpha=-slip;alpha<=slip;alpha++)
                   for (int beta=-slip;beta<=slip;beta++)
                        if(((i+alpha)>=1)&((j+beta)<=L)&((i+alpha)<=L)&((j+beta)>=1))
                            System.out.println("P "+Integer.toString(i+alpha)+" "+Integer.toString(j+beta)+" 1");
        }
        }
        catch (Exception e){
              System.err.println("Error: plot File Problem" + e.getMessage());
           }
    }
}

/*               for (int alpha=1;alpha<=slip;alpha++){
                   if((i-alpha)>=1){
               System.out.println("P "+Integer.toString(i-alpha)+" "+Integer.toString(j)+" 1");
                   }
               System.out.println("P "+Integer.toString(i+alpha)+" "+Integer.toString(j)+" 1");
               System.out.println("P "+Integer.toString(i)+" "+Integer.toString(j-alpha)+" 1");
               if((j+alpha)<=L){
               System.out.println("P "+Integer.toString(i)+" "+Integer.toString(j+alpha)+" 1");
*/