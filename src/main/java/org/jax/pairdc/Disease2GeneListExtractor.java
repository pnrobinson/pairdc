package org.jax.pairdc;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multimap;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.io.obo.hpo.HpoDiseaseAnnotationParser;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;

import java.util.Collection;
import java.util.List;
import java.util.Map;

public class Disease2GeneListExtractor {
    /** This is a version of the disease map that has OMIM diseases that correspond to the dated phenotype.hpoa file.*/
    private final Map<TermId, HpoDisease> diseaseMap;
    private final Multimap<TermId,TermId> geneToDiseaseMap;
    private final Map<TermId,String> geneIdToSymbolMap;

    public Disease2GeneListExtractor(String phenotypeDotHpoa,
                                     Ontology hpo,
                                     Multimap<TermId,TermId> gene2disease,
                                     Map<TermId,String> geneId2sym) {
        List<String> databases = ImmutableList.of("OMIM"); // restrict ourselves to OMIM entries
        this.diseaseMap = HpoDiseaseAnnotationParser.loadDiseaseMap(phenotypeDotHpoa, hpo,databases);
        geneToDiseaseMap = gene2disease;
        geneIdToSymbolMap = geneId2sym;
        init();
    }




    private void init() {
        int i=0;
        int j = 0;
        for (TermId geneId : geneToDiseaseMap.keys()) {
            Collection<TermId> diseaseIds = geneToDiseaseMap.get(geneId);
            for (TermId diseaseId : diseaseIds) {
                if (! diseaseMap.containsKey(diseaseId)) {
                   // HpoDisease d = diseaseMap.get(diseaseId);
                   // String name = d.getName();
                    String sym = geneIdToSymbolMap.get(geneId);
                    i++;
                    System.out.println(i + ") Could not find disease id for gene " + sym + " (" + diseaseId.getValue() +")" );
                } else {
                    j++;
                }
            }
        }
        System.out.println("Found " + j);


    }








}
