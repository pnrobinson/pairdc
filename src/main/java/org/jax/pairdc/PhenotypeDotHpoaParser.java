package org.jax.pairdc;

import java.io.*;
import java.text.ParseException;
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

        createDatedFile(date1, phenotypeDotHpoaFileOld);
        createDatedFile(date2, phenotypeDotHpoaFileNew);

    }


    /**
     * Output a version of the HPOA file that corresponds to date d, i.e., has no lines that
     * were created after d. Also, filter out non-OMIM lines
     * @param d target date
     * @param pathToNewHpoaFile path to phenotype.hpoa file
     */
    private void createDatedFile(Date d, String pathToNewHpoaFile) {
        int n_skipped_lines = 0;
        int n_accepted_lines = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(this.phenotypeDotHpoaFile));
            String line;
            BufferedWriter writer = new BufferedWriter(new FileWriter(pathToNewHpoaFile));
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.startsWith("DatabaseID")) {
                    writer.write(line + "\n"); // output header
                    continue;
                }
                String [] A = line.split("\t");
                if (A.length != 12) {
                    // should never happen
                    throw new RuntimeException("Malformed line " + line);
                }
                String datestring = A[11]; // can be single or multiple date.
                // we can just take the first date and not worry about modification dates
                int i = datestring.indexOf(";");
                if (i>0) {
                    datestring = datestring.substring(0, i);
                }
                i = datestring.indexOf("[");
                if (i>0) {
                    datestring = datestring.substring(i+1);
                }
                i = datestring.indexOf("]");
                if (i>0) {
                    datestring = datestring.substring(0,i);
                }
                Date biocurationDate = new SimpleDateFormat("yyyy-MM-DD").parse(datestring);
                if (biocurationDate.before(d)) {
                    writer.write(line + "\n");
                    n_accepted_lines ++;
                } else {
                    // otherwise the biocuration was after d, so we skip
                    n_skipped_lines++;
                }
            }
        } catch (IOException | ParseException e) {
            e.printStackTrace();
        }
        System.out.printf("Skipped %d lines that were curated after the target date of %s (wrote %d lines)\n",
                n_skipped_lines, d.toString(), n_accepted_lines);

    }



}
