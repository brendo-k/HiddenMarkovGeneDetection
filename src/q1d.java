import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Stack;

public class q1d {


    public static void main(String[] args){


        ArrayList<String[]> prediction = q1.getAnnotation("Question1c.txt");
        ArrayList<String[]> annotation = q1.getAnnotation("Vibrio_vulnificus.ASM74310v1.37.gff3");
        Stack<String[]> preStack = new Stack<>();
        Stack<String[]> annoStack = new Stack<>();

        preStack.addAll(prediction);
        annoStack.addAll(annotation);
        int same = 0;
        int start = 0;
        int end = 0;
        int sum = 0;


        while(!preStack.isEmpty() && !annoStack.isEmpty()){
            if(!preStack.peek()[0].equals(annoStack.peek()[0])){
                if(preStack.peek()[0].compareTo(annoStack.peek()[0]) > 0){
                    preStack.pop();
                }else{
                    annoStack.pop();
                }
            }else {
                String[] pred = preStack.peek();
                String[] anno = annoStack.peek();
                if (Integer.parseInt(pred[2]) <= Integer.parseInt(anno[1])) {
                    annoStack.pop();
                }else if(Integer.parseInt(anno[2]) <= Integer.parseInt(pred[1])){
                    preStack.pop();
                }else{
                    if(pred[1].equals(anno[1]) && pred[2].equals(anno[2])){
                        System.out.println(pred[1] + "\t" + pred[2]);
                        sum += Integer.parseInt(pred[2]) - Integer.parseInt(pred[1]);
                        same++;
                    }
                    else if(pred[1].equals(anno[1])){
                        start++;
                    }else if(pred[2].equals(2)){
                        end++;
                    }
                    preStack.pop();
                    annoStack.pop();
                }
            }
        }
        double len = (double)sum / same;
        System.out.println(len);
        System.out.println(same);
        System.out.println(start);
        System.out.println(end);

    }

}
