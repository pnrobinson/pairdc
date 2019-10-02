package org.jax.pairdc;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;


/** Parse the phenotype.hpoa file and create versions of the file
 * that correspond to date1 and date2.
 */
public class PhenotypeDotHpoaParser {


    private final Date date1;

    private final Date date2;

    private final File phenotypeDotHpoaFile;

    private String phenotypeDotHpoaFileOld;

    private String phenotypeDotHpoaFileNew;


    public PhenotypeDotHpoaParser(String path, Date d1, Date d2) {
        this.date1 = d1;
        this.date2 = d2;
        this.phenotypeDotHpoaFile = new File(path);
        if (! phenotypeDotHpoaFile.exists() ) {
            throw new RuntimeException("Could not find path to phenotype.hpoa file");
        }
    }


    public void createDatedPhenotypeHpoaFiles() {
        String dirname = phenotypeDotHpoaFile.getParent();
        String bname = phenotypeDotHpoaFile.getName();
        System.out.println("dir " + dirname + " base " + bname);
        String pattern = "yyyy-MM";
        SimpleDateFormat simpleDateFormat = new SimpleDateFormat(pattern);
        String basename1 = String.format("phenotype-%s.hpoa", simpleDateFormat.format(date1));
        phenotypeDotHpoaFileOld = String.format("%s%s%s", dirname,File.separator,basename1 );
        System.out.println(phenotypeDotHpoaFileOld);
        String basename2 = String.format("phenotype-%s.hpoa", simpleDateFormat.format(date2));
        phenotypeDotHpoaFileNew = String.format("%s%s%s", dirname,File.separator,basename2 );
        System.out.println(phenotypeDotHpoaFileNew);

    }



    private void createDatedFile(Date d, String basename) {

    }



}
