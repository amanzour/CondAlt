/*
 * This program reads a *.plot file and and a real number args[1]
 * It outputs a percentile (args[1]) of centroid to enforce (F) for hybrid-ss-min
 */
package selectcentroidbp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 *
 * @author manzouro
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        //float portion=Float.parseFloat(args[1]);
        int EXTEND=Integer.parseInt(args[2]);
        String seq="";
        try{
        BufferedReader im = new BufferedReader(new FileReader(args[0]));
        String str;
        while ((str = im.readLine()) != null)
            if (!str.startsWith(">"))
                seq=str;
        }        
        catch (Exception e){
              System.err.println("Error: Fasta File Problem" + e.getMessage());
           }    
        
        int counter=0;
        List<Integer> iList = new ArrayList<Integer>();
        List<Integer> jList = new ArrayList<Integer>();
        List<Float> pList = new ArrayList<Float>();
        try{
        BufferedReader in = new BufferedReader(new FileReader(args[1]));
        String str;
        int i=0;int j=0;float prob=0;
        float maxvalue=0;int iextra=0;int jextra=0;float probextra=0;
        while ((str = in.readLine()) != null){
            if ((!str.startsWith("i"))&&(Float.parseFloat(str.split("\t")[2])>0.5f)){
               i=Integer.parseInt(str.split("\t")[0]);
               j=Integer.parseInt(str.split("\t")[1]);
               prob=Float.parseFloat(str.split("\t")[2]);
               iList.add(i);jList.add(j);pList.add(prob);
               counter++;
        }
            if ((!str.startsWith("i"))&&(Float.parseFloat(str.split("\t")[2])>maxvalue)){
               iextra=Integer.parseInt(str.split("\t")[0]);
               jextra=Integer.parseInt(str.split("\t")[1]);
               probextra=Float.parseFloat(str.split("\t")[2]);
               maxvalue=Float.parseFloat(str.split("\t")[2]);
        }  
        }
        if (counter==0){
            iList.add(iextra);jList.add(jextra);pList.add(probextra);
            counter++;
        }
        }
        catch (Exception e){
              System.err.println("Error: plot File Problem" + e.getMessage());
           }
        
        int[] stem = new int[counter];
        int stemindex=1; stem[0]=stemindex;
        for (int cur=1;cur<counter;cur++)
            if (((iList.get(cur)-iList.get(cur-1)==1)&(jList.get(cur)-jList.get(cur-1)==-1))||((iList.get(cur)-iList.get(cur-1)==2)&(jList.get(cur)-jList.get(cur-1)==-1))||((iList.get(cur)-iList.get(cur-1)==1)&(jList.get(cur)-jList.get(cur-1)==-2)))
                stem[cur]=stemindex;
            else{
                stemindex++; 
                stem[cur]=stemindex;
            }
        int maxvalue=-1;
        int curmax=-1;
        int curval=0;
        int maxset=1;
      for (int i=1;i<=stemindex;i++){
          curval=0;
         for (int j=0;j<stem.length;j++){
            if (stem[j]==i)
                curval++;
         }
         if (curval>maxvalue){
            maxvalue=curval;
            maxset=i;
         }
      }
      int outerleft=-1;int outerright=-1;int initial=0;
      int innerleft=-1;int innerright=-1;
      for (int print=0;print<counter;print++)
          if (stem[print]==maxset){
              if (initial==0){
              outerleft=iList.get(print);outerright=jList.get(print);initial=1;
              }
             System.out.println("F "+iList.get(print)+" "+jList.get(print)+" 1"); 
             innerleft=iList.get(print);innerright=jList.get(print);
          }
      if (EXTEND>0){
        int counter2=0;
        List<Integer> iList2 = new ArrayList<Integer>();
        List<Integer> jList2 = new ArrayList<Integer>();
        List<Float> pList2 = new ArrayList<Float>();
      try{
        BufferedReader in2 = new BufferedReader(new FileReader(args[1]));
        String str;
        int i2=0;int j2=0;float prob2=0;
        
        while ((str = in2.readLine()) != null)
            if ((!str.startsWith("i"))&&(Float.parseFloat(str.split("\t")[2])>0.0f)){
               i2=Integer.parseInt(str.split("\t")[0]);
               j2=Integer.parseInt(str.split("\t")[1]);
               prob2=Float.parseFloat(str.split("\t")[2]);
               iList2.add(i2);jList2.add(j2);pList2.add(prob2);
               counter2++;
        }
        }
        catch (Exception e){
              System.err.println("Error: plot File Problem" + e.getMessage());
           }
           
  //extenting stem outwards 
      for (int print=0;print<counter2;print++){
          if ((iList2.get(print)==(outerleft-1))&&(jList2.get(print)==(outerright+1))){
              if (IsBP(seq.charAt(iList2.get(print)),seq.charAt(jList2.get(print))))
                   System.out.println("F "+iList2.get(print)+" "+jList2.get(print)+" 1");
         outerleft=iList2.get(print);outerright=jList2.get(print);print=0;
      }
      }
      
      //extending stem inwards     
         for (int print=0;print<counter2;print++){
         if ((iList2.get(print)==(innerleft+1))&&(jList2.get(print)==(innerright-1))){
             if (IsBP(seq.charAt(iList2.get(print)),seq.charAt(jList2.get(print))))
                System.out.println("F "+iList2.get(print)+" "+jList2.get(print)+" 1");
         innerleft=iList2.get(print);innerright=jList2.get(print);print=0;
      }
      }
    }
    }
    public static boolean IsBP(char i, char j){
        if (i=='A' & j=='U')
            return true;
        if (i=='U' & j=='A')
            return true;
        if (i=='C' & j=='G')
            return true;
        if (i=='G' & j=='C')
            return true;
        return false;
    }
}