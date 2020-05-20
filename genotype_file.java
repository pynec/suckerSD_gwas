package gwas_sim;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class genotype_file {
	int nloci;
	int nind;
	int nchr; 
	List<String> chr_list;
	List<String> chr_num;
	
	public genotype_file(int numind, int numloci, int numchr) {
		nind = numind;
		nloci = numloci;
		nchr = numchr;
	}
	public genotype_file() {}
	
	//method that creates the genotype file to be fed into GEMMA 
	//chr, minor, major, genotypes per individual (row = nind) 
	public List<List<String>> create_file(){
		List<List<String>> out_data = new ArrayList<>();
		//anno_list and chr_num_list = chr_list for the anno_file() method 
		List<String> anno_list = new ArrayList<>(); List<String> chr_num_list = new ArrayList<>();
		//list of possible genotypes for GEMMA 
		List<String> genotypes = new ArrayList<>();
		genotypes.add("0"); genotypes.add("1"); genotypes.add("2");
		String chr = "chr";
		//set up each row for the final output (chr, minor, major, genotypes) 
		for(int i = 0; i < nloci; i++) {
			List<String> nucleotides = new ArrayList<>();
			nucleotides.add("A"); nucleotides.add("T"); nucleotides.add("C"); nucleotides.add("G"); 
			List<String> out_data_row = new ArrayList<>();
			//choose a random number (from 1 to nchr) 
			Integer rand_char = (int)(Math.random()*((nchr-1)+1))+1;
			chr_num_list.add(rand_char.toString());
			// "chr" + random number = chr information (ex. chr3) 
			String chr_info = chr + rand_char.toString();	
			anno_list.add(chr_info);
			String minor = null; String major = null;
			//assign a random nucleotide for major and minor allele 
			for(int choose_nuc = 0; choose_nuc < 2; choose_nuc++) {
				Random rand = new Random();
				int index = rand.nextInt(nucleotides.size());
				if(choose_nuc == 0) { minor = nucleotides.get(index); }
				if(choose_nuc == 1) { major = nucleotides.get(index);}
				nucleotides.remove(index);
			}
			out_data_row.add(chr_info); out_data_row.add(minor); out_data_row.add(major);
			//get a list of random genotypes (one per individual) 
			genotype_file genotype_obj = new genotype_file();
			List<String> ind_geno = genotype_obj.get_genotypes(genotypes, nind);
			for(int j = 0; j < nind; j++) {
				out_data_row.add(ind_geno.get(j));
			}
			//add each row to the final output list 
			out_data.add(out_data_row);
		}
		chr_list = anno_list; chr_num = chr_num_list;
		return out_data;
	}
	
	//method that creates the phenotype file to be fed into GEMMA 
	//row = nind, col = number of phenotypes (sex = 2) 
	public List<Integer> phenotype_file(){
		List<Integer> out = new ArrayList<>();
		List<Integer> trait_val = new ArrayList<>();
		//sex = binary trait (0, 1);
		trait_val.add(0); trait_val.add(1);
		for(int i = 0; i < nind; i++) {
			Random r = new Random();
			int rIndex = r.nextInt(trait_val.size());
			out.add(trait_val.get(rIndex));
		}
		return out;	
	}
	
	//method that creates the annotation file to be fed into GEMMA 
	//row = nloci, col = SNP ID, base pair position, chr num 
	public List<List<String>> anno_file(){
		List<List<String>> anno = new ArrayList<>();
		//add a SNP_ID to chr number to make each one unique 
		for(int i = 0; i < nloci; i++) {
			List<String> anno_row = new ArrayList<>();
			int SNP_ID = (int)(Math.random()*((10-1)+1))+1;
			anno_row.add(chr_list.get(i) + SNP_ID);
			Integer bp_pos = (int)(Math.random()*((10000000-1000000)+1000000))+1000000;;
			anno_row.add(bp_pos.toString());
			anno_row.add(chr_num.get(i));
			anno.add(anno_row);
		}
		return anno;
		
	}
	
	//method that returns a list of random genotypes (0,1,2) per individual 
	public List<String> get_genotypes(List<String> genotypes, int nind) {
		List<String> return_list = new ArrayList<>();
		Random r = new Random();
		for(int i = 0; i < nind; i++) {
			int rIndex = r.nextInt(genotypes.size());
			return_list.add(genotypes.get(rIndex));}
		return return_list;	
	}
	
	//method that writes the output to a file 
	public static void write_to_file(List<List<String>> output) throws IOException {
		FileWriter writer = new FileWriter("geno_nosig.txt");
		for(List<String> row : output) {
			writer.write(row + System.lineSeparator());
		}
		writer.close();	
	}
	//method that writes the phenotype output to a file 
	public static void write_pheno_file(List<Integer> output) throws IOException {
		FileWriter writer = new FileWriter("pheno_nosig.txt");
		for(int i = 0; i < output.size(); i++) {
			writer.write(output.get(i) + System.lineSeparator());
		}
		writer.close();
	}
	//method that writes the annotation output to a file 
	public static void write_anno_file(List<List<String>> output) throws IOException {
		FileWriter writer = new FileWriter("anno_nosig.txt");
		for(List<String> row : output) {
			writer.write(row + System.lineSeparator());
		}
		writer.close();
		
	}

	public static void main(String[] args) throws Exception{
		//numind, numloci, numchr 
		genotype_file output = new genotype_file(30, 100, 10);
		List<List<String>> out = output.create_file();
		write_to_file(out);
		List<Integer>phenotypes = output.phenotype_file();
		write_pheno_file(phenotypes);
		List<List<String>> anno = output.anno_file();
		write_anno_file(anno);
		
	}

 }
