package org.jax.pairdc;

import org.monarchinitiative.phenol.ontology.data.TermId;

import java.util.Objects;

public class Gene2DiseaseAssociation {

    private final TermId geneId;
    private final TermId diseaseId;

    Gene2DiseaseAssociation(TermId gene, TermId disease){
        geneId = gene;
        diseaseId = disease;
    }


    public TermId getGeneId() {
        return geneId;
    }

    public TermId getDiseaseId() {
        return diseaseId;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Gene2DiseaseAssociation that = (Gene2DiseaseAssociation) o;
        return Objects.equals(geneId, that.geneId) &&
                Objects.equals(diseaseId, that.diseaseId);
    }

    @Override
    public int hashCode() {
        return Objects.hash(geneId, diseaseId);
    }


    @Override
    public String toString() {
        return getGeneToDiseaseAssociation();
    }

    String getGeneToDiseaseAssociation() {
        return String.format("%s\tg2d\t%s",geneId.getValue(), diseaseId.getValue() );
    }

//    public String getDiseaseToGeneAssociation() {
//        return String.format("%s\tdisease_to_gene\t%s", diseaseId.getValue(),geneId.getValue() );
//    }



}
