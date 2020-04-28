import java.io.*;
import java.util.*;

public class Viterbi {

    public static void main(String[] args) {
        String fastaPath = args[0];
        String configPath = args[1];

        //sequnces in the fasta file
        HashMap<String, String> seqs = q1.getSequence(fastaPath);
        //avg lengths or genic and intergenic
        double genicLength = 0;
        double intergenicLength = 0;

        //Emission tables
        HashMap<String, Double> genicEmissions = new HashMap<>();
        HashMap<Character, Double> intergenicEmission = new HashMap<>();
        HashMap<String, Double> startEmission = new HashMap<>();
        HashMap<String, Double> endEmission = new HashMap<>();

        //loading all info from viterbi file
        try {
            File file = new File(configPath);
            new FileWriter("Output.txt", false).close();
            BufferedReader br = new BufferedReader(new FileReader(file));
            String str;
            while ((str = br.readLine()) != null) {
                String[] values = str.split("\t");
                if (values[0].equals("Average gene size:")) {
                    genicLength = Double.parseDouble(values[1]);
                } else if (values[0].equals("Average intergenic size:")) {
                    intergenicLength = Double.parseDouble(values[1]);
                } else if (values[0].equals("Intergenic Emissions:")) {
                    for (int i = 1; i < values.length; i++) {
                        double value = Double.parseDouble(values[i].split(" ")[1]);
                        if(value == 0){
                            value = 0.001;
                        }else if(value == 1){
                            value = 0.999;
                        }
                        intergenicEmission.put(values[i].charAt(0), value);
                    }
                } else if (values[0].equals("Genic Emissions:")) {
                    for (int i = 1; i < values.length; i++) {
                        String[] keyValue = values[i].split(" ");
                        double value = Double.parseDouble(keyValue[1]);
                        if(value == 0){
                            value = 0.001;
                        }else if(value == 1){
                            value = 0.999;
                        }
                        genicEmissions.put(keyValue[0], value);
                    }
                } else if (values[0].equals("Start Emission:")) {
                    for (int i = 1; i < values.length; i++) {
                        String[] keyValue = values[i].split(" ");
                        double value = Double.parseDouble(keyValue[1]);
                        if(value == 0){
                            value = 0.001;
                        }else if(value == 1){
                            value = 0.999;
                        }
                        startEmission.put(keyValue[0], value);
                    }
                } else if (values[0].equals("End Emissions:")) {
                    for (int i = 1; i < values.length; i++) {
                        String[] keyValue = values[i].split(" ");
                        double value = Double.parseDouble(keyValue[1]);
                        if(value == 0){
                            value = 0.00001;
                        }else if(value == 1){
                            value = 0.999;
                        }
                        endEmission.put(keyValue[0], value);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        //making the transition table
        double[][] transition = new double[4][4];
        transition[0][1] = 1 / intergenicLength;
        transition[0][0] = (intergenicLength - 1) / intergenicLength;
        transition[2][2] = ((genicLength / 3) - 1) / (genicLength / 3);
        transition[2][3] = 1 / (genicLength / 3);
        transition[1][2] = 1;
        transition[3][0] = 1;

        //sorting the keys in the sequences so they're ordered alphabetically
        Set<String> seqKeys = (Set)seqs.keySet();
        ArrayList<String> keys = new ArrayList<>(seqKeys);
        Collections.sort(keys);

        //running viterbi on each sequence
        for (String seq : keys) {

            //intialize
            Pointer[][] viterbiArray = new Pointer[seqs.get(seq).length()][4];
            viterbiArray[0][0] = new Pointer(-1, -1, 1);
            for (int i = 1; i < viterbiArray[0].length; i++) {
                viterbiArray[0][i] = new Pointer(-1, -1, -Double.MAX_VALUE);
            }

            //filling viterbi array at position i and at all the states.
            //Pointers refrence where the value came from.
            //codon values are stored at the start of the the codon.
            for (int i = 1; i < viterbiArray.length; i++) {
                double[] intergenicScore = bestIntergenic(seqs.get(seq), viterbiArray, intergenicEmission, transition, i);
                //if came from stop codon point to the stop codon cell
                if(intergenicScore[1] == 3){
                    viterbiArray[i][0] = new Pointer((int) intergenicScore[1], i - 3,  intergenicScore[0]);
                }else {
                    viterbiArray[i][0] = new Pointer((int) intergenicScore[1], i - 1, intergenicScore[0]);
                }
                //checkng to make sure if you reference last codon or find current codon
                if (i - 1 > 3 && i + 3 < viterbiArray.length) {
                    double[] genicScore = bestCodon(seqs.get(seq), viterbiArray, genicEmissions, transition, i, 2);
                    viterbiArray[i][2] = new Pointer((int) genicScore[1],i - 3,  genicScore[0]);

                    double[] endScore = bestCodon(seqs.get(seq), viterbiArray, endEmission, transition, i, 3);
                    viterbiArray[i][3] = new Pointer((int) endScore[1],i - 3,  endScore[0]);
                } else {
                    viterbiArray[i][2] = new Pointer(-1, -1,-Double.MAX_VALUE);
                    viterbiArray[i][3] = new Pointer(-1, -1,-Double.MAX_VALUE);
                }
                //checking if you can look for a codon in the sequence
                if (i + 3 < viterbiArray.length) {
                    double[] startScore = bestStart(seqs.get(seq), viterbiArray, startEmission, transition, i);
                    viterbiArray[i][1] = new Pointer((int) startScore[1],i - 1,  startScore[0]);
                } else {
                    viterbiArray[i][1] = new Pointer(-1, -1, -Double.MAX_VALUE);
                }

            }
            //finding the state path from the viterbi array
            String path = backtrack(viterbiArray);

            //printing the appropriate information to a file
            try {
                BufferedWriter writer = new BufferedWriter(new FileWriter("Output.txt", true));
                int counter = 1;
                for(int i = 0; i < path.length(); i++){
                    if(path.charAt(i) == 'S'){
                        int[] range = getRange(path, i, counter);
                        counter += ((range[1]-range[0])/3+1)*2;
                        writer.write(seq + "\tena\tCDS\t" + range[0] + "\t" + range[1] + "\t.\t+\t0\t+\n");
                        writer.write("###\n");
                    }
                }
                writer.close();

            }catch (IOException e){
                e.printStackTrace();
            }
        }
    }

    //gets the range of a gene
    public static int[] getRange(String path, int start, int offset){
        for(int i = start+1; i < path.length(); i++){
            if(path.charAt(i) == 'I'){
                int difference = i-start;
                int[] range = {start+offset, start+difference*3-1+offset};
                return range;
            }
        }
        return null;
    }

    //backtracks viterbi array
    public static String backtrack(Pointer[][] array){
        int end = array.length-1;
        Pointer cur = array[end][0];
        int state = 0;
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < array[end].length; i++){
           // System.out.println(array[end][i].value);
            if(cur.value < array[end][i].value){
                state = i;
                cur = array[end][i];

            }
        }
        findState(sb, state);
        while(cur.endState != -1){
            state = cur.endState;
            findState(sb, state);
            cur = array[cur.endEmission][cur.endState];
        }
        return sb.toString();

    }

    //finds the state based on an int
    public static void findState(StringBuilder sb, int state){
        if(state == 0){
            sb.insert(0, "I");
        }else if(state == 1){
            sb.insert(0, "S");
        }else if(state == 2){
            sb.insert(0, "G");
        }else if(state == 3){
            sb.insert(0, "E");
        }
    }

    //finds the max prob for codon (used for genic and end state). Output is max[0] = value, max[1] = last state
    public static double[] bestCodon(String seq, Pointer[][] viterbiArray, HashMap<String, Double> codonEmission, double[][] transition, int i, int state){
        double[] possibleBest = new double[4];
        double[] max = new double[2];
        max[0] = -Double.MAX_VALUE;
        for(int j = 0; j < possibleBest.length; j++){
            String key = seq.substring(i, i+3);
            double emission;
            if(codonEmission.containsKey(key)) {
                emission = codonEmission.get(key);
            }else{
                emission = 0.00001;
            }
            Double prob;

            prob = viterbiArray[i-3][j].value + log(transition[j][state]) + log(emission);

            if(max[0] < prob){
                max[0] = prob;
                max[1] = j;
            }
        }
        return max;
    }

    //finds the max prob for intergenic. Output is max[0] = value, max[1] = last state
    public static double[] bestIntergenic(String seq, Pointer[][] viterbiArray, HashMap<Character, Double> intergenic, double[][] transition, int i){
        double[] possibleBest = new double[4];
        double[] max = new double[2];
        max[0] = -Double.MAX_VALUE;
        for(int j = 0; j < possibleBest.length; j++){
            Double prob;
            double weight = intergenic.get(seq.charAt(i));
            if(weight == 0){
                weight = 0.00001;
            }
            if(j == 3 && i-3 >= 0){
                prob = viterbiArray[i-3][j].value + log(transition[j][0]) + log(weight);
            }else {
                prob = viterbiArray[i - 1][j].value + log(transition[j][0]) + log(weight);
            }
            if(max[0] < prob){
                max[0] = prob;
                max[1] = j;
            }
        }
        return max;
    }

    //returns the log value. -inf = -max double
    public static double log(double number){
        double ans = Math.log(number);
        if(Double.isInfinite(ans)){
            ans = -Double.MAX_VALUE;
        }
        return ans;
    }

    //finds the max prob for start state. Output is max[0] = value, max[1] = last state
    public static double[] bestStart(String seq, Pointer[][] viterbiArray, HashMap<String, Double> codonEmission, double[][] transition, int i){
        double[] possibleBest = new double[4];
        double[] max = new double[2];
        max[0] = -Double.MAX_VALUE;
        for(int j = 0; j < possibleBest.length; j++){
            String key = seq.substring(i, i+3);
            double emission = 0;
            if(codonEmission.containsKey(key)) {
                emission = codonEmission.get(key);
            }else{
                emission = 0.00001;
            }
            Double prob = viterbiArray[i-1][j].value + log(transition[j][1]) + log(emission);
            if(max[0] < prob){
                max[0] = prob;
                max[1] = j;
            }else{

            }
        }
        return max;
    }
}

//helper class to backgrack viterbi array
class Pointer{
    double value;
    int endState;
    int endEmission;

    public Pointer(int endState, int endEmission, double value){
        this.endState = endState;
        this.endEmission = endEmission;
        this.value = value;

    }
}