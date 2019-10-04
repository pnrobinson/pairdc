package org.jax.pairdc;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multimap;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.io.obo.hpo.HpoDiseaseAnnotationParser;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;

import java.util.*;

public class Disease2GeneListExtractor {
    /** This is a version of the disease map that has OMIM diseases that correspond to the dated phenotype.hpoa file.*/
    private final Map<TermId, HpoDisease> diseaseMap;
    private final Multimap<TermId,TermId> geneToDiseaseMap;
    private final Map<TermId,String> geneIdToSymbolMap;

    private final Set<Gene2DiseaseAssociation> datedGeneToDiseaseSet;

    public Disease2GeneListExtractor(String phenotypeDotHpoa,
                                     Ontology hpo,
                                     Multimap<TermId,TermId> gene2disease,
                                     Map<TermId,String> geneId2sym) {
        List<String> databases = ImmutableList.of("OMIM"); // restrict ourselves to OMIM entries
        this.diseaseMap = HpoDiseaseAnnotationParser.loadDiseaseMap(phenotypeDotHpoa, hpo,databases);
        geneToDiseaseMap = gene2disease;
        geneIdToSymbolMap = geneId2sym;
        datedGeneToDiseaseSet = new HashSet<>();
        init();
    }

    Set<Gene2DiseaseAssociation> getDatedGeneToDiseaseSet() {
        return datedGeneToDiseaseSet;
    }

    private void init() {
        int i=0;
        int j = 0;
        for (TermId geneId : geneToDiseaseMap.keys()) {
            Collection<TermId> diseaseIds = geneToDiseaseMap.get(geneId);
            Set<TermId> alreadySeen = new HashSet<>(); // to avoid duplication
            for (TermId diseaseId : diseaseIds) {
                if (alreadySeen.contains(diseaseId)) {
                    continue;
                } else {
                    alreadySeen.add(diseaseId);
                }
                if ( diseaseMap.containsKey(diseaseId)) {
                    Gene2DiseaseAssociation g2d = new Gene2DiseaseAssociation(geneId,diseaseId);
                   datedGeneToDiseaseSet.add(g2d);
                } else {
                    j++;
                }
            }
        }
        System.out.println("Novel genes: n = " + j);
        System.out.println("Gene to disease association s known at date: n = " + datedGeneToDiseaseSet.size());

    }








}
