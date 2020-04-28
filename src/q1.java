import java.util.*;
import java.io.*;

public class q1 {
    //arg[0] = gene anotation file path
    //arg[1] = fasta sequence file path
    public static void main(String[] args){
        //returns the genes from CDS type in annotations
        ArrayList<String[]> listOfGene = getAnnotation(args[0]);
        //gets the sequence from fasta files
        HashMap<String, String> sequence = getSequence(args[1]);

        //calculuating the sum of gene sizes
        ArrayList<Integer> size = new ArrayList<>();
        for(String[] str: listOfGene){
            size.add(Integer.parseInt(str[2]) - Integer.parseInt(str[1])+1);
        }
        double sum = 0;
        for(int i = 0; i < size.size(); i++){
            sum += size.get(i);
        }
        sum /= (double)size.size();
        //printing the average size
        System.out.println("Average gene size:"+ "\t" + sum );

        //gets the sequence of genes and the sequence of intergenics
        ArrayList<String> geneSequence = new ArrayList<>();
        ArrayList<String> intergenicSequence = new ArrayList<>();
        String cur = listOfGene.get(0)[0];
        while(!listOfGene.isEmpty()) {
            int counter = 1;
            ArrayList<String> genic = new ArrayList<>();
            ArrayList<String> intergenic = new ArrayList<>();
            String[] currentGene = listOfGene.get(0);
            listOfGene.remove(0);
            while (counter < sequence.get(cur).length()) {
                cur = currentGene[0];
                String fullSeq = sequence.get(cur);

                while (counter < Integer.parseInt(currentGene[1])) {
                    intergenic.add(Character.toString(fullSeq.charAt(counter - 1)));
                    counter++;
                }
                intergenicSequence.add(listToString(intergenic));
                intergenic.clear();
                if(counter >= sequence.get(cur).length()){
                    break;
                }

                while (counter <= Integer.parseInt(currentGene[2])) {
                    genic.add(Character.toString(fullSeq.charAt(counter - 1)));
                    counter++;
                }
                geneSequence.add(listToString(genic));
                genic.clear();

                if (listOfGene.isEmpty() || !listOfGene.get(0)[0].equals(cur)) {
                    currentGene[1] = Integer.toString(sequence.get(cur).length());
                } else {
                    currentGene = listOfGene.get(0);
                    listOfGene.remove(0);
                }
            }
        }

        //summing the length of intergenic regions
        int intergenicSum = 0;
        for(String inter: intergenicSequence){
            intergenicSum += inter.length();
        }
        intergenicSum/= intergenicSequence.size();
        System.out.println("Average intergenic length: " + intergenicSum);

        //gets the emission occurences from genic and intergenic sequences
        HashMap<Character, Double> intergenicOccurence = intergenicOccurence(intergenicSequence);
        HashMap<String, Double> genicOccurence = genicOccurence(geneSequence);
        HashMap<String, Double> startOccurence = startOccurence(geneSequence);
        HashMap<String, Double> endOccurence  = endOccurence(geneSequence);

        //printing values
        System.out.println("Intergenic Occurrence:");
        for(Character car: intergenicOccurence.keySet()){
            System.out.print(car + " " + intergenicOccurence.get(car) + "  //");
        }
        System.out.println();
        System.out.println("Genic Occurrence:");
        for(String car: genicOccurence.keySet()){
            System.out.print(car + " " + genicOccurence.get(car) + "  //");
        }

        //writing important values for Viterbi.java into viterbi_param.txt file
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("viterbi_param.txt"));
            writer.write("Average gene size:"+ "\t" + sum + "\n");
            writer.write("Average intergenic size:"+ "\t" + intergenicSum + "\n");
            String str = "";
            for(Character car: intergenicOccurence.keySet()){
                str += (car + " " + intergenicOccurence.get(car) + "\t");
            }
            writer.write("Intergenic Emissions:\t" + str + "\n");
            str = "";
            for(String car: genicOccurence.keySet()){
                str += (car + " " + genicOccurence.get(car) + "\t");
            }
            writer.write("Genic Emissions:" + "\t" + str + "\n");
            str = "";
            for(String car: startOccurence.keySet()){
                str += (car + " " + startOccurence.get(car) + "\t");
            }
            writer.write("Start Emission:" + "\t" + str + "\n");
            str = "";
            for(String car: endOccurence.keySet()){
                str += (car + " " + endOccurence.get(car) + "\t");
            }
            writer.write("End Emissions:" + "\t" + str + "\n");

            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    //returns the lines CDC type genes in an arraylist
    public static ArrayList<String[]> getAnnotation(String fileName) {
        File file = new File(fileName);
        ArrayList<String[]> listOfGene = new ArrayList<>();
        try{
            BufferedReader br = new BufferedReader(new FileReader(file));

            String st = null;
            while((st = br.readLine()) != null){
                listOfGene.add(st.split("\t"));
            }
            Iterator<String[]> iterator = listOfGene.iterator();
            while(iterator.hasNext()){
                String[] stringarray = iterator.next();
                if(stringarray.length != 9){
                    iterator.remove();
                }else if(stringarray[6].equals("-")){
                    iterator.remove();
                }else if(!stringarray[2].equals("CDS")){
                    iterator.remove();
                }
            }
            for(int i = 0; i < listOfGene.size(); i++){
                listOfGene.set(i, Arrays.copyOfRange(listOfGene.get(i), 0, 5));
            }

            for(int i = 0; i  < listOfGene.size(); i++){
                String[] gene = listOfGene.get(i);
                String[] important = {gene[0], gene[3], gene[4]};
                listOfGene.set(i, important);
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return listOfGene;
    }

    //returns the dictionary of sequences.
    //key = sequence name, value = sequence.
    public static HashMap<String, String> getSequence(String fileName) {
        File file = new File(fileName);
        HashMap<String, String> sequence = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));

            String st = null;
            ArrayList<String> fullSequence = new ArrayList<>();
            st = br.readLine();
            String key = st.split(" ")[0].substring(1);
            while ((st = br.readLine()) != null) {
                if(st.charAt(0) == '>') {
                    sequence.put(key, listToString(fullSequence));
                    key = st.split(" ")[0].substring(1);
                    fullSequence.clear();
                }else{
                    fullSequence.add(st);
                }
            }
            sequence.put(key, q1.listToString(fullSequence));
            br.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
        return sequence;
    }

    //turns a list of strings to one string
    public static String listToString(ArrayList<String> str) {
        StringBuilder sb = new StringBuilder();

        // Appends characters one by one
        for (String ch : str) {
            sb.append(ch);
        }

        // convert in string
        String string = sb.toString();
        return string;
    }

    //gets the occurrences of intergenic nucleotides
    public static HashMap<Character, Double> intergenicOccurence(ArrayList<String> intergenic){
        HashMap<Character, Double> occurrence = new HashMap<Character, Double>() {{
            put('A', 0.0);
            put('C', 0.0);
            put('T', 0.0);
            put('G', 0.0);
        }};
            for(String seq: intergenic){
            for(Character nt: seq.toCharArray()){
                occurrence.replace(nt, occurrence.get(nt) + 1);
            }
        }
        int sum = 0;
        for(Double value: occurrence.values()){
            sum += value;
        }
        for(Character value: occurrence.keySet()){
            occurrence.replace(value, occurrence.get(value)/sum);
        }
        return occurrence;
    }

    //gets the occurrence of genic codons.
    public static HashMap<String, Double> genicOccurence(ArrayList<String> genics){
        HashMap<String, Double> occurrence = new HashMap<>();

        for(String seq: genics){
            int offset = seq.length()%3;
            for(int i = 0+offset; i < seq.length(); i+=3){
                String codon = "" + seq.charAt(i) + seq.charAt(i+1) + seq.charAt(i+2);
                if(!occurrence.containsKey(codon)){
                    occurrence.put(codon, 1.0);
                }else{
                    occurrence.replace(codon, occurrence.get(codon) + 1);
                }
            }
        }
        int sum = 0;
        for(Double value: occurrence.values()){
            sum += value;
        }
        for(String value: occurrence.keySet()){
            occurrence.replace(value, occurrence.get(value)/sum);
        }
        return occurrence;
    }

    //gets the occurrence of start codons
    public static HashMap<String, Double> startOccurence(ArrayList<String> genics){
        HashMap<String, Double> occurrence = new HashMap<>();

        for(String seq: genics){
            String codon = "" + seq.charAt(0) + seq.charAt(1) + seq.charAt(2);
            if(!occurrence.containsKey(codon)){
                occurrence.put(codon, 1.0);
            }else{
                occurrence.replace(codon, occurrence.get(codon) + 1);
            }
        }
        int sum = 0;
        for(Double value: occurrence.values()){
            sum += value;
        }
        for(String value: occurrence.keySet()){
            occurrence.replace(value, occurrence.get(value)/sum);
        }
        return occurrence;
    }

    //gets the occurence of stop codons
    public static HashMap<String, Double> endOccurence(ArrayList<String> genics){
        HashMap<String, Double> occurrence = new HashMap<>();

        for(String seq: genics){
            String codon = "" + seq.charAt(seq.length()-3) + seq.charAt(seq.length()-2) + seq.charAt(seq.length()-1);
            if(!occurrence.containsKey(codon)){
                occurrence.put(codon, 1.0);
            }else{
                occurrence.replace(codon, occurrence.get(codon) + 1);
            }
        }
        int sum = 0;
        for(Double value: occurrence.values()){
            sum += value;
        }
        for(String value: occurrence.keySet()){
            occurrence.replace(value, occurrence.get(value)/sum);
        }
        return occurrence;
    }


}
