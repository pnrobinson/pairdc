package org.jax.pairdc;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.io.OntologyLoader;
import org.monarchinitiative.phenol.io.assoc.HpoAssociationParser;
import org.monarchinitiative.phenol.io.obo.hpo.HpoDiseaseAnnotationParser;
import org.monarchinitiative.phenol.ontology.algo.InformationContentComputation;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;
import org.monarchinitiative.phenol.ontology.data.TermIds;
import org.monarchinitiative.phenol.ontology.similarity.PrecomputingPairwiseResnikSimilarity;
import org.monarchinitiative.phenol.ontology.similarity.ResnikSimilarity;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;


/**
 * To run this app, enter (adjusting the paths accordingly)
 * <pre>
 *   java -jar pairdc.jar  --hpo hp.obo -a phenotype.hpoa -o TEST --geneinfo Homo_sapiens_gene_info.gz --mim2genemedgen mim2gene_medgen
 * </pre>
 */
public class Main {
    /** Path to {@code hp.obo}. */
    @Parameter(names = {"-h","--hpo"}, description = "path to hp.obo file", required = true)
    private String hpoPath;
    /** Path to {@code phenotype.hpoa}. */
    @Parameter(names="-a", description = "path to phenotype.hpoa file", required = true)
    private String phenotypeDotHpoaPath;
    @Parameter(names="-o",description = "output file name")
    private String outname="pairwise_disease_similarity.tsv";
    /** Path to {@code Homo_sapiens_gene_info.gz} file. */
    @Parameter(names="--geneinfo",description = "path to downloaded file ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz")
    private String geneInfoPath;
    /** Path to {@code mim2gene_medgen} file with gene to disease associations.*/
    @Parameter(names="--mim2genemedgen",description = "path to downloaded file from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/mim2gene_medgen")
    private String mim2genMedgenPath;
    /** Number of threads to use. */
    private final int numThreads = 4;
    @Parameter(names="--date1", description = "First data for simulation, e.g., 4/2014", required = true)
    private String date1;
    @Parameter(names="--date2", description = "Second data for simulation, e.g., 9/2017", required = true)
    private String date2;
    @Parameter(names={"-t","--theshold"}, description = "SD threshold for disease similarity")
    private double diseaseSimilarityThreshold = 3.5;
    /** The gene to disease associations for date 1 (the "old data"). */
    private Set<Gene2DiseaseAssociation> dateOneg2dassocs;
    /** The gene to disease associations for date 2 (the "new data"). */
    private Set<Gene2DiseaseAssociation> dateTwog2dassocs;

    /** If true, perform pairwise gene-gene similarity analysis. Otherwise, perform pairwise disease-disease analysis.*/
    private boolean doGeneBasedAnalysis;

    private Ontology hpo;
    private Map<TermId, HpoDisease> diseaseMap;
    /** order list of HpoDiseases taken from the above map. */
    private List<HpoDisease> diseaseList;
    private Map<TermId,Integer> diseaseIdToIndexMap;
    private ResnikSimilarity resnikSimilarity;
    private  double[][] similarityScores;
    private Multimap<TermId,TermId> geneToDiseaseMap;
    private Map<TermId,String> geneIdToSymbolMap;
    private int n_diseases;

    private Main() { // no-op
    }

    static public void main(String[] args) {
        Main m = new Main();
        try {
            JCommander.newBuilder()
                    .addObject(m)
                    .build().
                    parse(args);
        } catch (ParameterException e) {
            e.printStackTrace();
        }
        try {
            m.run();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    /**
     * A gene may be associated with multiple diseases. Here, we take the maximum disease-disease
     * similarity between any disease associated with gene i and any disease associated with gene j
     * @param geneI the first disease gene
     * @param geneJ the second disease gene
     * @return maximum similarity between the genes
     */
    private double getMaximumGeneGeneSimilarity(TermId geneI, TermId geneJ) {
        double max = 0.0;
        Collection<TermId> diseasesI = this.geneToDiseaseMap.get(geneI);
        Collection<TermId> diseasesJ = this.geneToDiseaseMap.get(geneJ);
        if (diseasesI==null) {
            System.out.println("{Could not get diseases for gene " + geneI.getValue());
            return 0;
        }
        if (diseasesJ==null) {
            System.out.println("{Could not get diseases for gene " + geneJ.getValue());
            return 0;
        }
        for (TermId i : diseasesI) {
            for (TermId j : diseasesJ) {
                Integer index_i = this.diseaseIdToIndexMap.get(i);
                Integer index_j = this.diseaseIdToIndexMap.get(j);
                if (index_i==null) {
                    //System.err.println("[ERROR] COuld not retrieve index for disease " + i.getValue());
                    continue;
                }
                if (index_j==null) {
                    // System.err.println("[ERROR] Could not retrieve index for disease " + j.getValue());
                    continue;
                }
                double s = this.similarityScores[index_i][index_j];
                if (s>max) max=s;
            }
        }
        return max;
    }

    /**
     * Do an analysis to get the maximum pairwise similarity between genes, calculated on the basis
     * of phenotypic similarity of the diseases to which the genes are annotated.
     */
    private void performGeneBasedAnalysis() {
        HpoAssociationParser hpoAssociationParser = new HpoAssociationParser(this.geneInfoPath,this.mim2genMedgenPath,this.hpo);
        this.geneToDiseaseMap = hpoAssociationParser.getGeneToDiseaseIdMap();
        System.out.println("[INFO] geneToDiseaseMap with " + geneToDiseaseMap.size() + " entries");
        this.geneIdToSymbolMap = hpoAssociationParser.getGeneIdToSymbolMap();
        System.out.println("[INFO] geneIdToSymbolMap with " + geneIdToSymbolMap.size() + " entries");
        List<TermId> geneList = new ArrayList<>(geneToDiseaseMap.keySet());
        int N = geneList.size();
        double[][] geneSimilarityMatrix = new double[N][N];
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i=0;i<N-1;i++) {
            for (int j = i; j < N; j++) {
                try {
                    TermId geneI = geneList.get(i);
                    TermId geneJ = geneList.get(j);
                    if (geneI == null) {
                        System.err.println("gene i was null=" + i);
                        continue;
                    }
                    if (geneJ == null) {
                        System.err.println("gene j was null=" + j);
                        continue;
                    }
                    if (i > n_diseases) {
                        System.err.println(String.format("i=%d but n_diseases=%d", i, n_diseases));
                        continue;
                    }
                    if (j > n_diseases) {
                        System.err.println(String.format("j=%d but n_diseases=%d", j, n_diseases));
                        continue;
                    }
                    double sim = getMaximumGeneGeneSimilarity(geneI, geneJ);
                    geneSimilarityMatrix[i][j] = sim;
                    geneSimilarityMatrix[j][i] = sim;
                    stats.addValue(sim);
                } catch (Exception e){
                    System.err.println("i="+i+", j="+j+ " "+e.getMessage());
                }
            }
        }
        double mean = stats.getMean();
        double sd = stats.getStandardDeviation();
        System.out.println("\n\n[INFO] Done calculating gene based similarity matrix. Mean="+mean+", sd="+sd);
        double threshold = mean + 2.0*sd;
        System.out.println("[INFO] Writing pairwise gene similarity to file." );
        int aboveThreshold=0;
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(this.outname))){
            String [] fields = {"gene1","symbol1","gene2","symbol2","similarity"};
            String header = String.join("\t",fields);
            writer.write(header + "\n");
            for (int i=0;i<N-1;i++) {
                for (int j = i+1; j < N; j++) {
                    if (geneSimilarityMatrix[i][j]>threshold) {
                        TermId geneId1=geneList.get(i);
                        TermId geneId2=geneList.get(j);
                        String g1 = geneId1.getValue();
                        String g2 = geneList.get(j).getValue();
                        String symbol1 = this.geneIdToSymbolMap.get(geneId1);
                        String symbol2 = this.geneIdToSymbolMap.get(geneId2);

                        writer.write(g1 + "\t" + symbol1 + "\t" + g2+"\t"+symbol2+"\t"+geneSimilarityMatrix[i][j] + "\n");
                        aboveThreshold++;
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(String.format("[INFO] skipped vales: %d, good values %d",-42,stats.getN()));
        System.out.println(String.format("[INFO] Wrote %d above threshold (%.3f) pairwise interactions.",aboveThreshold,threshold) );
    }


    /**
     * Calculate the pairwise disease-disease similarities.
     * @param stats descriptive statices about the comparisons
     */
    private void performDiseaseBasedAnalysis(DescriptiveStatistics stats, Writer writer) throws  IOException {
        double mean = stats.getMean();
        double sd = stats.getStandardDeviation();
        System.out.println("\n\nMean="+mean+", sd="+sd);

        double threshold = mean + diseaseSimilarityThreshold*sd;
        System.out.printf("[INFO] Setting disease/disease similarity threshold to mean + %f*SD = %.2f\n", diseaseSimilarityThreshold,threshold);
        System.out.println("[INFO] Writing disease/disease pairwise phenotype similarity to file: \"" + outname +"\"" );
        int aboveThreshold=0;
            String [] fields = {"disease1","d2d","disease2","similarity"};
            String header = String.join("\t",fields);
            writer.write(header + "\n");
            int N = diseaseList.size();
            for (int i=0;i<N-1;i++) {
                for (int j = i+1; j < N; j++) {
                    if (similarityScores[i][j]>threshold) {
                        String d1 = diseaseList.get(i).getDiseaseDatabaseId().getValue();
                        String d2 = diseaseList.get(j).getDiseaseDatabaseId().getValue();
                        writer.write(d1 + "\td2d\t" + d2 + "\t" + similarityScores[i][j] + "\n");
                        aboveThreshold++;
                    }
                }
            }

        System.out.println(String.format("[INFO] Wrote %d above threshold (%.3f) pairwise interactions.",aboveThreshold,threshold) );
    }






    private void init() {
        HpoAssociationParser hpoAssociationParser = new HpoAssociationParser(this.geneInfoPath,this.mim2genMedgenPath,this.hpo);
        this.geneToDiseaseMap = hpoAssociationParser.getGeneToDiseaseIdMap();
        System.out.println("[INFO] geneToDiseaseMap with " + geneToDiseaseMap.size() + " entries");
        this.geneIdToSymbolMap = hpoAssociationParser.getGeneIdToSymbolMap();
        System.out.println("[INFO] geneIdToSymbolMap with " + geneIdToSymbolMap.size() + " entries");
        if (hpoPath==null || phenotypeDotHpoaPath == null) {
            System.err.println("[ERROR] Must pass path-to-hp.obo and path-to-phenotype.hpoa");
            System.exit(1);
        }
        this.hpo = OntologyLoader.loadOntology(new File(hpoPath));
        System.out.println("[INFO] DONE: Loading HPO");

        List<String> databases = ImmutableList.of("OMIM"); // restrict ourselves to OMIM entries
        this.diseaseMap = HpoDiseaseAnnotationParser.loadDiseaseMap(this.phenotypeDotHpoaPath, hpo,databases);
        System.out.println("[INFO] DONE: Loading phenotype.hpoa");

    }



    private void outputNewDiseaseGeneAssociations() {
        try {
            String dt = date2.replace("/","_");
            String outname =String.format("g2d_associations_test_%s.tsv",dt );
            BufferedWriter writer = new BufferedWriter(new FileWriter(outname));
            for (Gene2DiseaseAssociation g2d : dateTwog2dassocs) {
                if (!dateOneg2dassocs.contains(g2d)) {
                    writer.write(g2d.getGeneToDiseaseAssociation() + "\n");
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    /**
     * Run application.
     */
    private void run() throws  Exception {
        init();

        Date dt1 = new SimpleDateFormat("MM/yyyy").parse(date1);
        Date dt2 = new SimpleDateFormat("MM/yyyy").parse(date2);

        PhenotypeDotHpoaParser phparser = new PhenotypeDotHpoaParser(this.phenotypeDotHpoaPath,dt1,dt2);
        phparser.createDatedPhenotypeHpoaFiles();
        String oldPhenotypePath = phparser.getOldPhenotypeDotHpoaFile();
        String newPhenotypePath = phparser.getNewPhenotypeDotHpoaFile();
        System.out.println("Phenotype on " + dt1.toString());
        Disease2GeneListExtractor extractor = new Disease2GeneListExtractor(oldPhenotypePath,hpo,geneToDiseaseMap,geneIdToSymbolMap);
        dateOneg2dassocs = extractor.getDatedGeneToDiseaseSet();
        System.out.println("Phenotype on " + dt2.toString());
        Disease2GeneListExtractor extractor2 = new Disease2GeneListExtractor(newPhenotypePath,hpo,geneToDiseaseMap,geneIdToSymbolMap);
        dateTwog2dassocs = extractor2.getDatedGeneToDiseaseSet();

        System.out.println("Sanity check -- is everything from date 1 in date 2?");
        for (Gene2DiseaseAssociation g2d : dateOneg2dassocs) {
            if (! dateTwog2dassocs.contains(g2d)) {
                // should never happen, the following is an insanity check
                throw new RuntimeException("\"New\" annotation file did not contain " + g2d.toString() + " Some thing is wrong!");
            }
        }

        outputNewDiseaseGeneAssociations();
        // Now output disease-disease similarity for the new dataset.
        List<String> databases = ImmutableList.of("OMIM"); // restrict ourselves to OMIM entries
        Map<TermId, HpoDisease> dateTwoDiseaseMap = HpoDiseaseAnnotationParser.loadDiseaseMap(newPhenotypePath, hpo,databases);

        // Compute list of annoations and mapping from OMIM ID to term IDs.
        final Map<TermId, Collection<TermId>> diseaseIdToTermIds = new HashMap<>();
        final Map<TermId, Collection<TermId>> termIdToDiseaseIds = new HashMap<>();

        for (TermId diseaseId : dateTwoDiseaseMap.keySet()) {
            HpoDisease disease = dateTwoDiseaseMap.get(diseaseId);

            List<TermId> hpoTerms = disease.getPhenotypicAbnormalityTermIdList();
            diseaseIdToTermIds.putIfAbsent(diseaseId, new HashSet<>());
            // add term ancestors
            final Set<TermId> inclAncestorTermIds = TermIds.augmentWithAncestors(hpo, Sets.newHashSet(hpoTerms), true);
            for (TermId tid : inclAncestorTermIds) {
                termIdToDiseaseIds.putIfAbsent(tid, new HashSet<>());
                termIdToDiseaseIds.get(tid).add(diseaseId);
                diseaseIdToTermIds.get(diseaseId).add(tid);
            }
        }

        // Compute information content of HPO terms, given the term-to-disease annotation.
        System.out.println("[INFO] Performing IC precomputation...");
        final Map<TermId, Double> icMap =
                new InformationContentComputation(hpo)
                        .computeInformationContent(termIdToDiseaseIds);
        System.out.println("[INFO] DONE: Performing IC precomputation");


        // Initialize Resnik similarity precomputation
        System.out.println("[INFO] Performing Resnik precomputation...");
        final PrecomputingPairwiseResnikSimilarity pairwiseResnikSimilarity =
                new PrecomputingPairwiseResnikSimilarity(hpo, icMap, numThreads);
        System.out.println("[INFO] DONE: Performing Resnik precomputation");
        this.resnikSimilarity =
                new ResnikSimilarity(pairwiseResnikSimilarity, false);
        System.out.println(String.format("name: %s  params %s",
                resnikSimilarity.getName(),
                resnikSimilarity.getParameters()));
        System.out.println("[INFO] Calculating pairwise phenotype similarity for " + diseaseMap.size() + " diseases." );

        this.diseaseList = new ArrayList<>(diseaseMap.values());
        this.diseaseIdToIndexMap=new HashMap<>();
        for (int i=0;i<diseaseList.size();i++){
            this.diseaseIdToIndexMap.put(diseaseList.get(i).getDiseaseDatabaseId(),i);
        }
        n_diseases = diseaseList.size();
        int expectedTotal = n_diseases*(n_diseases-1)/2;
        this.similarityScores = new double[n_diseases][n_diseases];
        DescriptiveStatistics stats = new DescriptiveStatistics();
        int c=0;

        for (int i=0;i<n_diseases;i++) {
            for (int j=i;j<n_diseases;j++) {
                HpoDisease d1 = diseaseList.get(i);
                HpoDisease d2 = diseaseList.get(j);
                List<TermId> pheno1 = d1.getPhenotypicAbnormalityTermIdList();
                List<TermId> pheno2 = d2.getPhenotypicAbnormalityTermIdList();
                double similarity = resnikSimilarity.computeScore(pheno1, pheno2);
                if (Double.isNaN(similarity)) similarity = 0;
                similarityScores[i][j]=similarity;
                similarityScores[j][i]=similarity; // symmetric
                stats.addValue(similarity);
                if (++c%10000==0) {
                    System.out.print(String.format("Got %d/%d similarity counts (%.1f%%)\r",c,expectedTotal,100.0*(double)c/expectedTotal));
                }
            }
        }

        System.out.println(String.format("[INFO] Disease analysis: n similarity values %d", stats.getN()));
        if (doGeneBasedAnalysis) {
            performGeneBasedAnalysis();
        } else {
            try {
                BufferedWriter writer = new BufferedWriter(new FileWriter(this.outname));
                performDiseaseBasedAnalysis(stats, writer);
                for (Gene2DiseaseAssociation g2d : dateOneg2dassocs ){
                    writer.write(g2d.getGeneToDiseaseAssociation() + "\n");
                }
            } catch (IOException e) {
            e.printStackTrace();
        }
        }
    }

}
