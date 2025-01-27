package evolution.speciation;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import sr.evolution.tree.SRTree;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;
import sr.evolution.sranges.StratigraphicRange;
import ssr.speciation.SRangesBirthDeathSkylineModel;


import java.util.ArrayList;

/**
 * Created by gavryusa on 24/07/17.
 * Edited by Kate Truman, January 2025
 */
public class SRangesBirthDeathSkylineModelTest extends TestCase {

    @Test
    public void testLikelihood() throws Exception {
//String newick = "(((t8_4:0.2304798054,t8_3:0):0.3794035359,t8_2:0):0.6130463995,(t1_1:0.05680109053,((t2_6:0.3317049952,((((t5_10:0.06002887006,t9_1:0.06213831523):0.1608192763,t3_1:0.6675525472):0.1112413445,((((t6_16:0,t6_15:0):0.009338884954,(t7_18:0,t7_17:0):0.009338884954):0.3025040342,(t12_12:0,t12_11:0):0.3118429191):0.3369682326,((t10_1:0.1907707396,(t11_14:0,t11_13:0):0.1907707396):0.4123371679,(t4_8:0,t4_7:0):0.6031079075):0.04570324422):0.12998274):0.6050714199,t5_9:0):0.3889511814):0.7453878975,t2_5:0):0.3771968665):1.024216084):0.08038265895;";

        String newick = "(((((A:3.4,2_last:0.0):1.0,2_first:0.0):0.7,(B:3.5,(3_last:1.7,3_first:0.0):0.8):1.6):0.55,1_last:0.0):0.85,1_first:0.0):0.5";
        Tree tree_initial = new TreeParser(newick, false);
        StratigraphicRange sr1 = new StratigraphicRange();
        Taxon taxon1_first = new Taxon("1_first");
        Taxon taxon1_last = new Taxon("1_last");
        sr1.setInputValue("firstOccurrence", taxon1_first);
        sr1.setInputValue("lastOccurrence", taxon1_last);
        StratigraphicRange sr2 = new StratigraphicRange();
        Taxon taxon2_first = new Taxon("2_first");
        Taxon taxon2_last = new Taxon("2_last");
        sr2.setInputValue("firstOccurrence", taxon2_first);
        sr2.setInputValue("lastOccurrence", taxon2_last);
        StratigraphicRange sr3 = new StratigraphicRange();
        Taxon taxon3_first = new Taxon("3_first");
        Taxon taxon3_last = new Taxon("3_last");
        sr3.setInputValue("firstOccurrence", taxon3_first);
        sr3.setInputValue("lastOccurrence", taxon3_last);
        ArrayList<StratigraphicRange> sranges = new ArrayList<>();
        sranges.add(sr1);
        sranges.add(sr2);
        sranges.add(sr3);
        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.assignFrom(tree_initial);

        SRangesBirthDeathSkylineModel model = new SRangesBirthDeathSkylineModel();
        model.setInputValue("tree", tree);
        model.setInputValue("intervalTimes", new RealParameter("6 5 4 3 2 1 0"));
//        model.setInputValue("reverseTimeArrays", new BooleanParameter("true true true true true"));
        model.setInputValue("origin", new RealParameter("7.0"));
//        model.setInputValue("birthRate", new RealParameter("1.5"));
//        model.setInputValue("deathRate", new RealParameter("0.5"));
//        model.setInputValue("samplingRate", new RealParameter("0.1"));
//        model.setInputValue("removalProbability", new RealParameter("0.0"));
        model.setInputValue("removalProbability", new RealParameter("0.0 0.0 0.0 0.0 0.0 0.0 0.0"));
        model.setInputValue("rho", new RealParameter("0.5"));

        model.setInputValue("netDiversification", new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"));
        model.setInputValue("turnOver", new RealParameter("0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333"));
        model.setInputValue("samplingProportion", new RealParameter("0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667"));
        model.setInputValue("contemp", true);
        model.initAndValidate();
        assertEquals(-33.74668640318646, model.calculateLogP(), 1e-7);


    }

    @Test
    public void testLikelihoodSingleTip() throws Exception {
        String newick = "(((((A:3.4,2_last:0.0):1.0,2_first:0.0):0.7,(B:3.5,(3_first:0.0):0.8):1.6):0.55,1_last:0.0):0.85,1_first:0.0):0.5";

        Tree tree_initial = new TreeParser(newick, false);
        StratigraphicRange sr1 = new StratigraphicRange();
        Taxon taxon1_first = new Taxon("1_first");
        Taxon taxon1_last = new Taxon("1_last");
        sr1.setInputValue("firstOccurrence", taxon1_first);
        sr1.setInputValue("lastOccurrence", taxon1_last);
        StratigraphicRange sr2 = new StratigraphicRange();
        Taxon taxon2_first = new Taxon("2_first");
        Taxon taxon2_last = new Taxon("2_last");
        sr2.setInputValue("firstOccurrence", taxon2_first);
        sr2.setInputValue("lastOccurrence", taxon2_last);
        StratigraphicRange sr3 = new StratigraphicRange();
        Taxon taxon3_first = new Taxon("3_first");
        sr3.setInputValue("firstOccurrence", taxon3_first);
        sr3.setInputValue("lastOccurrence", taxon3_first);
        ArrayList<StratigraphicRange> sranges = new ArrayList<>();
        sranges.add(sr1);
        sranges.add(sr2);
        sranges.add(sr3);
        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.assignFrom(tree_initial);

        SRangesBirthDeathSkylineModel model = new SRangesBirthDeathSkylineModel();
        model.setInputValue("tree", tree);
        model.setInputValue("intervalTimes", new RealParameter("6 5 4 3 2 1 0"));
        model.setInputValue("origin", new RealParameter("7.0"));
        model.setInputValue("removalProbability", new RealParameter("0.0 0.0 0.0 0.0 0.0 0.0 0.0"));
        model.setInputValue("rho", new RealParameter("0.5"));
        model.setInputValue("netDiversification", new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"));
        model.setInputValue("turnOver", new RealParameter("0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333"));
        model.setInputValue("samplingProportion", new RealParameter("0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667"));
        model.setInputValue("contemp", true);
        model.initAndValidate();
        assertEquals(-30.588945157591503, model.calculateLogP(), 1e-7);

    }

    @Test
    public void testLikelihoodSkylineRates() throws Exception {
        String newick = "(((((A:3.4,2_last:0.0):1.0,2_first:0.0):0.7,(B:3.5,(3_last:1.7,3_first:0.0):0.8):1.6):0.55,1_last:0.0):0.85,1_first:0.0):0.5";

        Tree tree_initial = new TreeParser(newick, false);
        StratigraphicRange sr1 = new StratigraphicRange();
        Taxon taxon1_first = new Taxon("1_first");
        Taxon taxon1_last = new Taxon("1_last");
        sr1.setInputValue("firstOccurrence", taxon1_first);
        sr1.setInputValue("lastOccurrence", taxon1_last);
        StratigraphicRange sr2 = new StratigraphicRange();
        Taxon taxon2_first = new Taxon("2_first");
        Taxon taxon2_last = new Taxon("2_last");
        sr2.setInputValue("firstOccurrence", taxon2_first);
        sr2.setInputValue("lastOccurrence", taxon2_last);
        StratigraphicRange sr3 = new StratigraphicRange();
        Taxon taxon3_first = new Taxon("3_first");
        Taxon taxon3_last = new Taxon("3_last");
        sr3.setInputValue("firstOccurrence", taxon3_first);
        sr3.setInputValue("lastOccurrence", taxon3_last);
        ArrayList<StratigraphicRange> sranges = new ArrayList<>();
        sranges.add(sr1);
        sranges.add(sr2);
        sranges.add(sr3);
        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.assignFrom(tree_initial);

        SRangesBirthDeathSkylineModel model = new SRangesBirthDeathSkylineModel();
        model.setInputValue("tree", tree);
        model.setInputValue("intervalTimes", new RealParameter("4.7 3.5 1 0"));
        model.setInputValue("origin", new RealParameter("7.0"));
        model.setInputValue("removalProbability", new RealParameter("0.0 0.0 0.0 0.0"));
        model.setInputValue("rho", new RealParameter("0.7"));
        model.setInputValue("netDiversification", new RealParameter("0.25 0.5 0.09 0.5"));
        model.setInputValue("turnOver", new RealParameter("0.1666666667 0.2857142857 0.1 0.375"));
        model.setInputValue("samplingProportion", new RealParameter("0.6666666667 0.6 0.9756097561 0.7692307692"));
        model.setInputValue("contemp", true);
        model.initAndValidate();
        assertEquals(-24.68765411973161, model.calculateLogP(), 1e-7);
    }
}
