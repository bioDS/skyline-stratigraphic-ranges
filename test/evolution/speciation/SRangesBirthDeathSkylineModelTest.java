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
        model.setInputValue("conditionOnSurvival", true);
        model.setInputValue("conditionOnRhoSampling", false);
        model.setInputValue("conditionOnRoot", false);



        model.initAndValidate();
        double ans = model.calculateLogP();
        System.out.println(ans);
        // Condition on rho sampling
        //assertEquals(-33.34091738060979, model.calculateLogP(), 1e-9);
        // Conditioning on survival
        assertEquals(-33.38391796963602, model.calculateLogP(), 1e-9);




    }

    @Test
    public void testLikelihoodSingletonSR() throws Exception {
//String newick = "(((t8_4:0.2304798054,t8_3:0):0.3794035359,t8_2:0):0.6130463995,(t1_1:0.05680109053,((t2_6:0.3317049952,((((t5_10:0.06002887006,t9_1:0.06213831523):0.1608192763,t3_1:0.6675525472):0.1112413445,((((t6_16:0,t6_15:0):0.009338884954,(t7_18:0,t7_17:0):0.009338884954):0.3025040342,(t12_12:0,t12_11:0):0.3118429191):0.3369682326,((t10_1:0.1907707396,(t11_14:0,t11_13:0):0.1907707396):0.4123371679,(t4_8:0,t4_7:0):0.6031079075):0.04570324422):0.12998274):0.6050714199,t5_9:0):0.3889511814):0.7453878975,t2_5:0):0.3771968665):1.024216084):0.08038265895;";

        String newick = "((((A:3.4,2_last:0.0):1.0,2_first:0.0):0.7,(B:3.5,(3_last:1.7,3_first:0.0):0.8):1.6):0.55,1_first:0.0):1.35";

        Tree tree_initial = new TreeParser(newick, false);
        StratigraphicRange sr1 = new StratigraphicRange();
        Taxon taxon1_first = new Taxon("1_first");
        sr1.setInputValue("firstOccurrence", taxon1_first);
        sr1.setInputValue("lastOccurrence", taxon1_first);
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
        model.setInputValue("origin", new RealParameter("7.0"));
        model.setInputValue("removalProbability", new RealParameter("0.0 0.0 0.0 0.0 0.0 0.0 0.0"));
        model.setInputValue("rho", new RealParameter("0.5"));

        model.setInputValue("netDiversification", new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"));
        model.setInputValue("turnOver", new RealParameter("0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333"));
        model.setInputValue("samplingProportion", new RealParameter("0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667"));
        model.setInputValue("contemp", true);
        model.setInputValue("conditionOnSurvival", false);
        model.setInputValue("conditionOnRhoSampling", false);
        model.setInputValue("conditionOnRoot", false);
        model.initAndValidate();
        double ans = model.calculateLogP();
        System.out.println(ans);
        assertEquals(-31.141006127536, ans, 1e-9);


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
        model.setInputValue("conditionOnSurvival", false);
        model.setInputValue("conditionOnRhoSampling", false);
        model.setInputValue("conditionOnRoot", false);
        model.initAndValidate();
        System.out.println(model.calculateLogP());
        assertEquals(-30.588945157591503, model.calculateLogP(), 1e-9);
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
        model.setInputValue("conditionOnSurvival", false);
        model.setInputValue("conditionOnRhoSampling", false);
        model.setInputValue("conditionOnRoot", false);
        model.initAndValidate();
        System.out.println(model.calculateLogP());
        assertEquals(-24.68765411973161, model.calculateLogP(), 1e-7);
    }

    @Test
    public void testLikelihoodSkylineRatesTwo() throws Exception {
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
        model.setInputValue("intervalTimes", new RealParameter("5.1 4.7 3.5 1 0"));
        model.setInputValue("origin", new RealParameter("7.0"));
        model.setInputValue("removalProbability", new RealParameter("0.0 0.0 0.0 0.0 0.0"));
        model.setInputValue("rho", new RealParameter("0.5"));
        model.setInputValue("netDiversification", new RealParameter("1.0 1.0 1.0 1.0 1.0"));
        model.setInputValue("turnOver", new RealParameter("0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333"));
        model.setInputValue("samplingProportion", new RealParameter("0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667"));
        model.setInputValue("contemp", true);

        model.setInputValue("conditionOnSurvival", false);
        model.setInputValue("conditionOnRhoSampling", false);
        model.setInputValue("conditionOnRoot", false);

        model.initAndValidate();
        double ans = model.calculateLogP();
        System.out.println(ans);
        assertEquals(-33.74668640318646,ans , 1e-7);
    }


    @Test
    public void testLikelihoodSkylineRatesBDS() throws Exception {
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
        model.setInputValue("intervalTimes", new RealParameter("5.1 4.7 3.5 1 0"));
        model.setInputValue("origin", new RealParameter("7.0"));
        model.setInputValue("removalProbability", new RealParameter("0.0 0.0 0.0 0.0 0.0"));
        model.setInputValue("rho", new RealParameter("0.5"));
        model.setInputValue("birthRate", new RealParameter("1.5 1.5 1.5 1.5 1.5"));
        model.setInputValue("deathRate", new RealParameter("0.5 0.5 0.5 0.5 0.5"));
        model.setInputValue("samplingRate", new RealParameter("0.1 0.1 0.1 0.1 0.1"));
        model.setInputValue("contemp", true);

        model.setInputValue("conditionOnSurvival", false);
        model.setInputValue("conditionOnRhoSampling", false);
        model.setInputValue("conditionOnRoot", false);

        model.initAndValidate();
        double ans = model.calculateLogP();
        System.out.println(ans);
        assertEquals(-33.74668640318646,ans , 1e-7);
    }


    @Test
    public void testLikelihoodLphy() throws Exception {
        String newick = "((f1_first:0.0,((f10_first:0.0,(f23_first:0.0,f25_first:0.08375786698920429):0.048842486639731675):0.3517911055607139,f1_last:0.002883633520253337):1.1422401632612886):0.0811674323290883,(f3_first:0.0,(f3_last:0.0,((f5_first:0.0,f5_last:0.12842802311274903):0.26008920211757003,(((((f14_first:0.0,(((((f42_first:0.0,(f42_last:0.0,((f62_first:0.0,(((f113_first:0.0,f113_last:0.049538460083955044):0.2695067740941405,f62_last:0.3190452341780955):5.699051691705481E-4,(f112_first:0.0,f112_last:0.22683808448199477):0.09277705486527127):0.2364943717326783):0.6049872532632139,(f89_first:0.0,f89_last:0.2361348002306598):0.17954812374800266):0.4235729034881264):0.14074985327708034):0.163637508160688,((((f68_first:0.0,f68_last:0.32604134291425124):0.21726117558596303,f53_first:0.11690355082080295):0.029605680696134895,(f66_first:0.0,((f114_first:0.3165583996670492,(f102_first:0.0,f102_last:0.27780939549871575):0.038749004168333434):0.26529989914647506,(f66_last:0.4146829297334672,(f107_first:0.0,f107_last:0.3526427252405628):0.06204020449290443):0.16717536908005703):0.3342358127819324):0.44343567446980514):0.05107197901205618,((f64_first:0.0,((f109_first:0.0697137210048998,(f121_first:0.0,f121_last:0.051595346811140724):0.01811837419375907):0.28866406939297695,(f64_last:0.06280557258969384,f122_first:0.06280557258969384):0.29557221780818294):0.5174029863746684):0.41326742174199893,(((f77_first:0.0,f77_last:0.7031702680034307):0.07047300546112345,(f101_first:0.0,(f101_last:0.3329710954311305,(f110_first:0.0,f110_last:0.1713896365425588):0.16158145888857173):0.23637777481152006):0.20429440322190362):0.23971853487924344,(f92_first:0.0,(f120_first:0.030545279819267943,f124_first:0.030545279819267943):0.7507430534115984):0.2320734751129313):0.2756863901707465):0.12155356656277383):0.478455264191735):0.20826813523222154,f31_first:0.24814782804777757):0.08268672664145438,(((((f30_first:0.0,((f104_first:0.49101842499555287,f30_last:0.49101842499555287):0.06494758598731915,(f103_first:0.2585548457845621,f115_first:0.2585548457845621):0.2974111651983099):0.012296738356025627):0.4811943947950563,(f91_first:0.0,f91_last:0.09916329697022275):0.11122076415769078):0.19186975530385442,(((f76_first:0.0,f76_last:0.10569766126405615):0.5816569917051787,f99_first:0.13268564729673715):0.35703826941090877,(((f98_first:0.0,f98_last:0.12188309026566303):0.5815885386058682,(f97_first:0.0,f97_last:0.3902693875267839):0.31320224134474733):0.013235947396091352,f90_first:0.3211167198191775):0.327685346112521):0.1969339770576648):0.19517305851279354,((((f74_first:0.0,f74_last:0.018242547203483728):0.12371937271195177,f86_first:0.13075037171694726):0.08794007462162523,(f79_first:0.0,(f108_first:0.404279590498846,f79_last:0.404279590498846):0.19308623505078282):0.5867599304626465):0.09912497404072074,f67_first:0.0011758831791903734):0.15324922789760587):0.5582501988700626,(f38_first:0.0,f38_last:0.033626886977095705):0.16668724667206258):0.1852617343220644):0.06866880236514605,(f33_first:0.0,f33_last:0.16829170274580418):0.2390182169225712):0.06419599899133566):0.21432142889729233,f7_first:0.028563327958678553):0.08005458521160325,((f12_first:0.0,f12_last:0.006006816316410646):0.050961278081560035,(f29_first:0.0,(((f46_first:0.0,((f60_first:0.0,((f93_first:0.0,(((f106_first:0.0,(f106_last:0.14322296504495946,f118_first:0.14322296504495946):0.185013064453447):0.1002822262817617,(f93_last:0.03292469547950541,f123_first:0.03292469547950541):0.3955935603006627):0.023803669501758273,(f105_first:0.174770675859228,f117_first:0.174770675859228):0.27755124942269843):0.45067911384591214):0.09075613528779514,f60_last:0.9937571744156337):0.14853972985434383):0.36634283355039154,f46_last:0.010181337999765061):0.06750856422943952):0.2859810762043795,f75_first:0.7132882835268204):0.0142956181919347,(f29_last:0.0,(((f69_first:0.0,f69_last:0.979704975149927):0.006060288955954896,f94_first:0.19507724851527974):0.3403383440556821,(f70_first:0.0,f70_last:0.7365739243412226):0.5895296838203414):0.5270649286793192):0.023256459605239588):0.25614216914187504):0.022646298946602528):0.45203924207350576):0.1406001216792454,(f20_first:0.0,f20_last:0.12372182123368969):0.3862181773150688):0.015292416863728509,(f16_first:0.0,(f16_last:0.1776500972375059,f21_first:0.11509821383363272):0.05044060122631322):0.3208233733716428):0.3258006617998168):0.13550825163158642):0.5824774868295002):0.07058604612793307):0.0";
        Tree tree_initial = new TreeParser(newick, false);
        StratigraphicRange sr0 = new StratigraphicRange();
        Taxon taxon0_first = new Taxon("f1_first");
        Taxon taxon0_last = new Taxon("f1_last");
        sr0.setInputValue("firstOccurrence", taxon0_first);
        sr0.setInputValue("lastOccurrence", taxon0_last);

        StratigraphicRange sr1 = new StratigraphicRange();
        Taxon taxon1_first = new Taxon("f10_first");
        sr1.setInputValue("firstOccurrence", taxon1_first);
        sr1.setInputValue("lastOccurrence", taxon1_first);

        StratigraphicRange sr2 = new StratigraphicRange();
        Taxon taxon2_first = new Taxon("f23_first");
        sr2.setInputValue("firstOccurrence", taxon2_first);
        sr2.setInputValue("lastOccurrence", taxon2_first);

        StratigraphicRange sr3 = new StratigraphicRange();
        Taxon taxon3_first = new Taxon("f25_first");
        sr3.setInputValue("firstOccurrence", taxon3_first);
        sr3.setInputValue("lastOccurrence", taxon3_first);

        StratigraphicRange sr4 = new StratigraphicRange();
        Taxon taxon4_first = new Taxon("f3_first");
        Taxon taxon4_last = new Taxon("f3_last");
        sr4.setInputValue("firstOccurrence", taxon4_first);
        sr4.setInputValue("lastOccurrence", taxon4_last);

        StratigraphicRange sr5 = new StratigraphicRange();
        Taxon taxon5_first = new Taxon("f5_first");
        Taxon taxon5_last = new Taxon("f5_last");
        sr5.setInputValue("firstOccurrence", taxon5_first);
        sr5.setInputValue("lastOccurrence", taxon5_last);

        StratigraphicRange sr6 = new StratigraphicRange();
        Taxon taxon6_first = new Taxon("f14_first");
        sr6.setInputValue("firstOccurrence", taxon6_first);
        sr6.setInputValue("lastOccurrence", taxon6_first);

        StratigraphicRange sr7 = new StratigraphicRange();
        Taxon taxon7_first = new Taxon("f42_first");
        Taxon taxon7_last = new Taxon("f42_last");
        sr7.setInputValue("firstOccurrence", taxon7_first);
        sr7.setInputValue("lastOccurrence", taxon7_last);

        StratigraphicRange sr8 = new StratigraphicRange();
        Taxon taxon8_first = new Taxon("f62_first");
        Taxon taxon8_last = new Taxon("f62_last");
        sr8.setInputValue("firstOccurrence", taxon8_first);
        sr8.setInputValue("lastOccurrence", taxon8_last);

        StratigraphicRange sr9 = new StratigraphicRange();
        Taxon taxon9_first = new Taxon("f113_first");
        Taxon taxon9_last = new Taxon("f113_last");
        sr9.setInputValue("firstOccurrence", taxon9_first);
        sr9.setInputValue("lastOccurrence", taxon9_last);

        StratigraphicRange sr10 = new StratigraphicRange();
        Taxon taxon10_first = new Taxon("f112_first");
        Taxon taxon10_last = new Taxon("f112_last");
        sr10.setInputValue("firstOccurrence", taxon10_first);
        sr10.setInputValue("lastOccurrence", taxon10_last);

        StratigraphicRange sr11 = new StratigraphicRange();
        Taxon taxon11_first = new Taxon("f89_first");
        Taxon taxon11_last = new Taxon("f89_last");
        sr11.setInputValue("firstOccurrence", taxon11_first);
        sr11.setInputValue("lastOccurrence", taxon11_last);

        StratigraphicRange sr12 = new StratigraphicRange();
        Taxon taxon12_first = new Taxon("f68_first");
        Taxon taxon12_last = new Taxon("f68_last");
        sr12.setInputValue("firstOccurrence", taxon12_first);
        sr12.setInputValue("lastOccurrence", taxon12_last);

        StratigraphicRange sr13 = new StratigraphicRange();
        Taxon taxon13_first = new Taxon("f53_first");
        sr13.setInputValue("firstOccurrence", taxon13_first);
        sr13.setInputValue("lastOccurrence", taxon13_first);

        StratigraphicRange sr14 = new StratigraphicRange();
        Taxon taxon14_first = new Taxon("f66_first");
        Taxon taxon14_last = new Taxon("f66_last");
        sr14.setInputValue("firstOccurrence", taxon14_first);
        sr14.setInputValue("lastOccurrence", taxon14_last);

        StratigraphicRange sr15 = new StratigraphicRange();
        Taxon taxon15_first = new Taxon("f114_first");
        sr15.setInputValue("firstOccurrence", taxon15_first);
        sr15.setInputValue("lastOccurrence", taxon15_first);

        StratigraphicRange sr16 = new StratigraphicRange();
        Taxon taxon16_first = new Taxon("f102_first");
        Taxon taxon16_last = new Taxon("f102_last");
        sr16.setInputValue("firstOccurrence", taxon16_first);
        sr16.setInputValue("lastOccurrence", taxon16_last);

        StratigraphicRange sr17 = new StratigraphicRange();
        Taxon taxon17_first = new Taxon("f107_first");
        Taxon taxon17_last = new Taxon("f107_last");
        sr17.setInputValue("firstOccurrence", taxon17_first);
        sr17.setInputValue("lastOccurrence", taxon17_last);

        StratigraphicRange sr18 = new StratigraphicRange();
        Taxon taxon18_first = new Taxon("f64_first");
        Taxon taxon18_last = new Taxon("f64_last");
        sr18.setInputValue("firstOccurrence", taxon18_first);
        sr18.setInputValue("lastOccurrence", taxon18_last);

        StratigraphicRange sr19 = new StratigraphicRange();
        Taxon taxon19_first = new Taxon("f109_first");
        sr19.setInputValue("firstOccurrence", taxon19_first);
        sr19.setInputValue("lastOccurrence", taxon19_first);

        StratigraphicRange sr20 = new StratigraphicRange();
        Taxon taxon20_first = new Taxon("f121_first");
        Taxon taxon20_last = new Taxon("f121_last");
        sr20.setInputValue("firstOccurrence", taxon20_first);
        sr20.setInputValue("lastOccurrence", taxon20_last);

        StratigraphicRange sr21 = new StratigraphicRange();
        Taxon taxon21_first = new Taxon("f122_first");
        sr21.setInputValue("firstOccurrence", taxon21_first);
        sr21.setInputValue("lastOccurrence", taxon21_first);

        StratigraphicRange sr22 = new StratigraphicRange();
        Taxon taxon22_first = new Taxon("f77_first");
        Taxon taxon22_last = new Taxon("f77_last");
        sr22.setInputValue("firstOccurrence", taxon22_first);
        sr22.setInputValue("lastOccurrence", taxon22_last);

        StratigraphicRange sr23 = new StratigraphicRange();
        Taxon taxon23_first = new Taxon("f101_first");
        Taxon taxon23_last = new Taxon("f101_last");
        sr23.setInputValue("firstOccurrence", taxon23_first);
        sr23.setInputValue("lastOccurrence", taxon23_last);

        StratigraphicRange sr24 = new StratigraphicRange();
        Taxon taxon24_first = new Taxon("f110_first");
        Taxon taxon24_last = new Taxon("f110_last");
        sr24.setInputValue("firstOccurrence", taxon24_first);
        sr24.setInputValue("lastOccurrence", taxon24_last);

        StratigraphicRange sr25 = new StratigraphicRange();
        Taxon taxon25_first = new Taxon("f92_first");
        sr25.setInputValue("firstOccurrence", taxon25_first);
        sr25.setInputValue("lastOccurrence", taxon25_first);

        StratigraphicRange sr26 = new StratigraphicRange();
        Taxon taxon26_first = new Taxon("f120_first");
        sr26.setInputValue("firstOccurrence", taxon26_first);
        sr26.setInputValue("lastOccurrence", taxon26_first);

        StratigraphicRange sr27 = new StratigraphicRange();
        Taxon taxon27_first = new Taxon("f124_first");
        sr27.setInputValue("firstOccurrence", taxon27_first);
        sr27.setInputValue("lastOccurrence", taxon27_first);

        StratigraphicRange sr28 = new StratigraphicRange();
        Taxon taxon28_first = new Taxon("f31_first");
        sr28.setInputValue("firstOccurrence", taxon28_first);
        sr28.setInputValue("lastOccurrence", taxon28_first);

        StratigraphicRange sr29 = new StratigraphicRange();
        Taxon taxon29_first = new Taxon("f30_first");
        Taxon taxon29_last = new Taxon("f30_last");
        sr29.setInputValue("firstOccurrence", taxon29_first);
        sr29.setInputValue("lastOccurrence", taxon29_last);

        StratigraphicRange sr30 = new StratigraphicRange();
        Taxon taxon30_first = new Taxon("f104_first");
        sr30.setInputValue("firstOccurrence", taxon30_first);
        sr30.setInputValue("lastOccurrence", taxon30_first);

        StratigraphicRange sr31 = new StratigraphicRange();
        Taxon taxon31_first = new Taxon("f103_first");
        sr31.setInputValue("firstOccurrence", taxon31_first);
        sr31.setInputValue("lastOccurrence", taxon31_first);

        StratigraphicRange sr32 = new StratigraphicRange();
        Taxon taxon32_first = new Taxon("f115_first");
        sr32.setInputValue("firstOccurrence", taxon32_first);
        sr32.setInputValue("lastOccurrence", taxon32_first);

        StratigraphicRange sr33 = new StratigraphicRange();
        Taxon taxon33_first = new Taxon("f91_first");
        Taxon taxon33_last = new Taxon("f91_last");
        sr33.setInputValue("firstOccurrence", taxon33_first);
        sr33.setInputValue("lastOccurrence", taxon33_last);

        StratigraphicRange sr34 = new StratigraphicRange();
        Taxon taxon34_first = new Taxon("f76_first");
        Taxon taxon34_last = new Taxon("f76_last");
        sr34.setInputValue("firstOccurrence", taxon34_first);
        sr34.setInputValue("lastOccurrence", taxon34_last);

        StratigraphicRange sr35 = new StratigraphicRange();
        Taxon taxon35_first = new Taxon("f99_first");
        sr35.setInputValue("firstOccurrence", taxon35_first);
        sr35.setInputValue("lastOccurrence", taxon35_first);

        StratigraphicRange sr36 = new StratigraphicRange();
        Taxon taxon36_first = new Taxon("f98_first");
        Taxon taxon36_last = new Taxon("f98_last");
        sr36.setInputValue("firstOccurrence", taxon36_first);
        sr36.setInputValue("lastOccurrence", taxon36_last);

        StratigraphicRange sr37 = new StratigraphicRange();
        Taxon taxon37_first = new Taxon("f97_first");
        Taxon taxon37_last = new Taxon("f97_last");
        sr37.setInputValue("firstOccurrence", taxon37_first);
        sr37.setInputValue("lastOccurrence", taxon37_last);

        StratigraphicRange sr38 = new StratigraphicRange();
        Taxon taxon38_first = new Taxon("f90_first");
        Taxon taxon38_last = new Taxon("f90_first");
        sr38.setInputValue("firstOccurrence", taxon38_first);
        sr38.setInputValue("lastOccurrence", taxon38_last);

        StratigraphicRange sr39 = new StratigraphicRange();
        Taxon taxon39_first = new Taxon("f74_first");
        Taxon taxon39_last = new Taxon("f74_last");
        sr39.setInputValue("firstOccurrence", taxon39_first);
        sr39.setInputValue("lastOccurrence", taxon39_last);

        StratigraphicRange sr40 = new StratigraphicRange();
        Taxon taxon40_first = new Taxon("f86_first");
        sr40.setInputValue("firstOccurrence", taxon40_first);
        sr40.setInputValue("lastOccurrence", taxon40_first);

        StratigraphicRange sr41 = new StratigraphicRange();
        Taxon taxon41_first = new Taxon("f79_first");
        Taxon taxon41_last = new Taxon("f79_last");
        sr41.setInputValue("firstOccurrence", taxon41_first);
        sr41.setInputValue("lastOccurrence", taxon41_last);

        StratigraphicRange sr42 = new StratigraphicRange();
        Taxon taxon42_first = new Taxon("f108_first");
        sr42.setInputValue("firstOccurrence", taxon42_first);
        sr42.setInputValue("lastOccurrence", taxon42_first);

        StratigraphicRange sr43 = new StratigraphicRange();
        Taxon taxon43_first = new Taxon("f67_first");
        sr43.setInputValue("firstOccurrence", taxon43_first);
        sr43.setInputValue("lastOccurrence", taxon43_first);

        StratigraphicRange sr44 = new StratigraphicRange();
        Taxon taxon44_first = new Taxon("f38_first");
        Taxon taxon44_last = new Taxon("f38_last");
        sr44.setInputValue("firstOccurrence", taxon44_first);
        sr44.setInputValue("lastOccurrence", taxon44_last);

        StratigraphicRange sr45 = new StratigraphicRange();
        Taxon taxon45_first = new Taxon("f33_first");
        Taxon taxon45_last = new Taxon("f33_last");
        sr45.setInputValue("firstOccurrence", taxon45_first);
        sr45.setInputValue("lastOccurrence", taxon45_last);

        StratigraphicRange sr46 = new StratigraphicRange();
        Taxon taxon46_first = new Taxon("f7_first");
        sr46.setInputValue("firstOccurrence", taxon46_first);
        sr46.setInputValue("lastOccurrence", taxon46_first);

        StratigraphicRange sr47 = new StratigraphicRange();
        Taxon taxon47_first = new Taxon("f12_first");
        Taxon taxon47_last = new Taxon("f12_last");
        sr47.setInputValue("firstOccurrence", taxon47_first);
        sr47.setInputValue("lastOccurrence", taxon47_last);

        StratigraphicRange sr48 = new StratigraphicRange();
        Taxon taxon48_first = new Taxon("f29_first");
        Taxon taxon48_last = new Taxon("f29_last");
        sr48.setInputValue("firstOccurrence", taxon48_first);
        sr48.setInputValue("lastOccurrence", taxon48_last);

        StratigraphicRange sr49 = new StratigraphicRange();
        Taxon taxon49_first = new Taxon("f46_first");
        Taxon taxon49_last = new Taxon("f46_last");
        sr49.setInputValue("firstOccurrence", taxon49_first);
        sr49.setInputValue("lastOccurrence", taxon49_last);

        StratigraphicRange sr50 = new StratigraphicRange();
        Taxon taxon50_first = new Taxon("f60_first");
        Taxon taxon50_last = new Taxon("f60_last");
        sr50.setInputValue("firstOccurrence", taxon50_first);
        sr50.setInputValue("lastOccurrence", taxon50_last);

        StratigraphicRange sr51 = new StratigraphicRange();
        Taxon taxon51_first = new Taxon("f93_first");
        Taxon taxon51_last = new Taxon("f93_last");
        sr51.setInputValue("firstOccurrence", taxon51_first);
        sr51.setInputValue("lastOccurrence", taxon51_last);

        StratigraphicRange sr52 = new StratigraphicRange();
        Taxon taxon52_first = new Taxon("f106_first");
        Taxon taxon52_last = new Taxon("f106_last");
        sr52.setInputValue("firstOccurrence", taxon52_first);
        sr52.setInputValue("lastOccurrence", taxon52_last);

        StratigraphicRange sr53 = new StratigraphicRange();
        Taxon taxon53_first = new Taxon("f118_first");
        sr53.setInputValue("firstOccurrence", taxon53_first);
        sr53.setInputValue("lastOccurrence", taxon53_first);

        StratigraphicRange sr54 = new StratigraphicRange();
        Taxon taxon54_first = new Taxon("f123_first");
        sr54.setInputValue("firstOccurrence", taxon54_first);
        sr54.setInputValue("lastOccurrence", taxon54_first);

        StratigraphicRange sr55 = new StratigraphicRange();
        Taxon taxon55_first = new Taxon("f105_first");
        sr55.setInputValue("firstOccurrence", taxon55_first);
        sr55.setInputValue("lastOccurrence", taxon55_first);

        StratigraphicRange sr56 = new StratigraphicRange();
        Taxon taxon56_first = new Taxon("f117_first");
        sr56.setInputValue("firstOccurrence", taxon56_first);
        sr56.setInputValue("lastOccurrence", taxon56_first);

        StratigraphicRange sr57 = new StratigraphicRange();
        Taxon taxon57_first = new Taxon("f75_first");
        sr57.setInputValue("firstOccurrence", taxon57_first);
        sr57.setInputValue("lastOccurrence", taxon57_first);

        StratigraphicRange sr58 = new StratigraphicRange();
        Taxon taxon58_first = new Taxon("f69_first");
        Taxon taxon58_last = new Taxon("f69_last");
        sr58.setInputValue("firstOccurrence", taxon58_first);
        sr58.setInputValue("lastOccurrence", taxon58_last);

        StratigraphicRange sr59 = new StratigraphicRange();
        Taxon taxon59_first = new Taxon("f94_first");
        sr59.setInputValue("firstOccurrence", taxon59_first);
        sr59.setInputValue("lastOccurrence", taxon59_first);

        StratigraphicRange sr60 = new StratigraphicRange();
        Taxon taxon60_first = new Taxon("f70_first");
        Taxon taxon60_last = new Taxon("f70_last");
        sr60.setInputValue("firstOccurrence", taxon60_first);
        sr60.setInputValue("lastOccurrence", taxon60_last);

        StratigraphicRange sr61 = new StratigraphicRange();
        Taxon taxon61_first = new Taxon("f20_first");
        Taxon taxon61_last = new Taxon("f20_last");
        sr61.setInputValue("firstOccurrence", taxon61_first);
        sr61.setInputValue("lastOccurrence", taxon61_last);

        StratigraphicRange sr62 = new StratigraphicRange();
        Taxon taxon62_first = new Taxon("f16_first");
        Taxon taxon62_last = new Taxon("f16_last");
        sr62.setInputValue("firstOccurrence", taxon62_first);
        sr62.setInputValue("lastOccurrence", taxon62_last);

        StratigraphicRange sr63 = new StratigraphicRange();
        Taxon taxon63_first = new Taxon("f21_first");
        sr63.setInputValue("firstOccurrence", taxon63_first);
        sr63.setInputValue("lastOccurrence", taxon63_first);

        ArrayList<StratigraphicRange> sranges = new ArrayList<>();
        sranges.add(sr0);
        sranges.add(sr1);
        sranges.add(sr2);
        sranges.add(sr3);
        sranges.add(sr4);
        sranges.add(sr5);
        sranges.add(sr6);
        sranges.add(sr7);
        sranges.add(sr8);
        sranges.add(sr9);
        sranges.add(sr10);
        sranges.add(sr11);
        sranges.add(sr12);
        sranges.add(sr13);
        sranges.add(sr14);
        sranges.add(sr15);
        sranges.add(sr16);
        sranges.add(sr17);
        sranges.add(sr18);
        sranges.add(sr19);
        sranges.add(sr20);
        sranges.add(sr21);
        sranges.add(sr22);
        sranges.add(sr23);
        sranges.add(sr24);
        sranges.add(sr25);
        sranges.add(sr26);
        sranges.add(sr27);
        sranges.add(sr28);
        sranges.add(sr29);
        sranges.add(sr30);
        sranges.add(sr31);
        sranges.add(sr32);
        sranges.add(sr33);
        sranges.add(sr34);
        sranges.add(sr35);
        sranges.add(sr36);
        sranges.add(sr37);
        sranges.add(sr38);
        sranges.add(sr39);
        sranges.add(sr40);
        sranges.add(sr41);
        sranges.add(sr42);
        sranges.add(sr43);
        sranges.add(sr44);
        sranges.add(sr45);
        sranges.add(sr46);
        sranges.add(sr47);
        sranges.add(sr48);
        sranges.add(sr49);
        sranges.add(sr50);
        sranges.add(sr51);
        sranges.add(sr52);
        sranges.add(sr53);
        sranges.add(sr54);
        sranges.add(sr55);
        sranges.add(sr56);
        sranges.add(sr57);
        sranges.add(sr58);
        sranges.add(sr59);
        sranges.add(sr60);
        sranges.add(sr61);
        sranges.add(sr62);
        sranges.add(sr63);
        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.assignFrom(tree_initial);

        SRangesBirthDeathSkylineModel model = new SRangesBirthDeathSkylineModel();
        model.setInputValue("tree", tree);
        model.setInputValue("intervalTimes", new RealParameter("3 2 1 0"));
        model.setInputValue("origin", new RealParameter("4.0"));
        model.setInputValue("removalProbability", new RealParameter("0.0 0.0 0.0 0.0"));
        model.setInputValue("rho", new RealParameter("0.9842431891153070"));
        model.setInputValue("netDiversification", new RealParameter("0.8796392547066560 0.8561234544571510 0.8650893306785030 0.8107079829462880 "));
        model.setInputValue("turnOver", new RealParameter("0.4196549226008420 0.7866014138753280 0.7258701963134100 0.3727589081350500"));
        model.setInputValue("samplingProportion", new RealParameter("0.6763761721058830 0.4517988068732030 0.25541481906884800 0.7789621706540440"));
        model.initAndValidate();
        double ans = model.calculateLogP();
        System.out.println(ans);


    }

    @Test
    public void testLikelihoodSLargeTree() throws Exception {
        String newick = "((t74_1:0.5183676227,((t38_1:0.06721284018,(t57_1:0.2091827123,t54_1:0.5444231388):0.2663953874):0.4169615277,((t67_1:0.06071098638,t44_1:0.07841508391):0.02535103175,(((((t8_1:0.02656962038,t46_1:0.1446107331):0.1664626977,((((t55_1:0.46663225,(t1_35:0.2498651272,t1_33:0):0.2167671227):0.1124660273,((t53_59:0.09887199938,t52_60:0.09887199938):0.1633361457,t53_57:0):0.3168901322):0.2988107614,t49_43:0):0.08664353278,((((t39_38:0.3995456205,t65_39:0.3995456205):0.07608750604,((t68_54:0.02247076183,t68_52:0):0.1482946044,t14_55:0.1707653662):0.3048677603):0.1067393297,((t23_45:0.2537365066,t16_47:0.2537365066):0.1066013398,t75_41:0.3603378465):0.2220346098):0.3034695115,t39_36:0):0.07871060374):0.5409178623):0.2445417434,t77_1:0.08390988151):0.03641923026,(t7_1:0.07939586374,t61_1:0.2879464697):0.2099953164):0.303817933,t62_5:2.09024934):0.86383299):0.08146975462):0.1769225896):0.5160259013,(t2_1:0.8107987027,(((t19_1:0.1709951962,(((((t30_27:0.5910243129,((t33_71:0.02026818522,t5_1:0.02026818522):0.1843990877,t33_69:0):0.3863570399):0.169037002,t45_22:0.447981802):0.872806652,(t43_1:0.005153761383,t31_1:0.2872419196):0.0182915175):0.1657281295,(t34_1:0.01319812915,t47_1:0.5644161765):0.4085346187):0.1845815716,t3_14:0):0.1842404518):0.7133656847,(((((t12_4:0.4198126475,((((t41_13:1.488931735,(((t37_65:0.05318776089,t9_67:0.05318776089):0.7870539546,(t35_29:0.5852354949,((((t66_64:0.06407216964,(t13_72:0.0005822999635,t27_74:0.0005822999635):0.06348986968):0.1098372975,t22_1:0.1726218266):0.1766046554,t66_62:0):0.1427396358,t50_31:0.4288205444):0.0919817366):0.2550062205):0.5799782334,(t64_1:0.1376593949,((t42_16:1.190803256,t21_1:0.006811195953):0.04670139866,((t20_1:0.1797930163,t60_1:0.4128313573):0.169122223,t28_20:0.8065558372):0.4309488171):0.003824879357):0.1788904151):0.06871178604):0.07128650743,(((((t48_26:0.7494946056,t40_1:0.01070018054):0.01054057277,t48_24:0):0.1592450521,t70_18:0.9192802305):0.2483028663,t51_1:0.1293601267):0.246278286,((t32_51:0.1746055724,t59_1:0.1746055724):0.1895445951,t32_49:0):1.049711215):0.1463568595):0.07016040752,t41_11:0):0.5002920884,t76_1:0.6534471572):0.2644277844):0.1155693866,(t73_1:0.5750588836,((t6_1:0.1684733718,t11_9:0):0.4034568713,(t15_1:0.9793441046,(t63_1:0.1104186289,t26_7:0.02564097301):0.1985442919):0.1692909344):0.03477073447):0.1240307332):0.1969817617,t12_2:0):0.02901065049,t10_1:0.1233843662):0.02291368474,(t71_1:0.2512826659,t36_1:0.1681021358):0.1402059386):0.1212097985):0.4497600315,(t58_1:0.3309153432,(t4_1:0.3268538318,(((t25_1:0.2168007838,((t72_1:0.1749694054,t18_1:0.3611771581):0.1642162904,(t29_1:0.02328003506,t24_1:0.825280175):0.09833783613):0.1906261948):0.05139049964,(t17_1:0.9791325493,t69_1:0.1173910335):0.2025095811):0.003110596624,t56_1:0.09679021826):0.1012050304):0.8779442821):0.008916652653):0.3274540253):0.07050271465):0.2714994241;";

        Tree tree_initial = new TreeParser(newick, false);
        StratigraphicRange sr1 = new StratigraphicRange();
        Taxon taxon1_first = new Taxon("t39_36");
        Taxon taxon1_last = new Taxon("t39_38");
        sr1.setInputValue("firstOccurrence", taxon1_first);
        sr1.setInputValue("lastOccurrence", taxon1_first);
        StratigraphicRange sr2 = new StratigraphicRange();
        Taxon taxon2_first = new Taxon("t53_57");
        Taxon taxon2_last = new Taxon("t53_59");
        sr2.setInputValue("firstOccurrence", taxon2_first);
        sr2.setInputValue("lastOccurrence", taxon2_last);
        StratigraphicRange sr3 = new StratigraphicRange();
        Taxon taxon3_first = new Taxon("t68_52");
        Taxon taxon3_last = new Taxon("t68_54");
        sr3.setInputValue("firstOccurrence", taxon3_first);
        sr3.setInputValue("lastOccurrence", taxon3_last);
        StratigraphicRange sr4 = new StratigraphicRange();
        Taxon taxon4_first = new Taxon("t33_69");
        Taxon taxon4_last = new Taxon("t33_71");
        sr4.setInputValue("firstOccurrence", taxon4_first);
        sr4.setInputValue("lastOccurrence", taxon4_last);
        StratigraphicRange sr5 = new StratigraphicRange();
        Taxon taxon5_first = new Taxon("t3_14");
        sr5.setInputValue("firstOccurrence", taxon5_first);
        sr5.setInputValue("lastOccurrence", taxon5_first);
        StratigraphicRange sr6 = new StratigraphicRange();
        Taxon taxon6_first = new Taxon("t1_33");
        Taxon taxon6_last = new Taxon("t1_35");
        sr6.setInputValue("firstOccurrence", taxon6_first);
        sr6.setInputValue("lastOccurrence", taxon6_last);
        StratigraphicRange sr7 = new StratigraphicRange();
        Taxon taxon7_first = new Taxon("t49_43");
        sr7.setInputValue("firstOccurrence", taxon7_first);
        sr7.setInputValue("lastOccurrence", taxon7_first);
        StratigraphicRange sr8 = new StratigraphicRange();
        Taxon taxon8_first = new Taxon("t11_9");
        sr8.setInputValue("firstOccurrence", taxon8_first);
        sr8.setInputValue("lastOccurrence", taxon8_first);
        StratigraphicRange sr9 = new StratigraphicRange();
        Taxon taxon9_first = new Taxon("t12_2");
        Taxon taxon9_last = new Taxon("t12_4");
        sr9.setInputValue("firstOccurrence", taxon9_first);
        sr9.setInputValue("lastOccurrence", taxon9_last);
        StratigraphicRange sr10 = new StratigraphicRange();
        Taxon taxon10_first = new Taxon("t41_11");
        Taxon taxon10_last = new Taxon("t41_13");
        sr10.setInputValue("firstOccurrence", taxon10_first);
        sr10.setInputValue("lastOccurrence", taxon10_last);
        StratigraphicRange sr11 = new StratigraphicRange();
        Taxon taxon11_first = new Taxon("t48_24");
        Taxon taxon11_last = new Taxon("t48_26");
        sr11.setInputValue("firstOccurrence", taxon11_first);
        sr11.setInputValue("lastOccurrence", taxon11_last);
        StratigraphicRange sr12 = new StratigraphicRange();
        Taxon taxon12_first = new Taxon("t66_62");
        Taxon taxon12_last = new Taxon("t66_64");
        sr12.setInputValue("firstOccurrence", taxon12_first);
        sr12.setInputValue("lastOccurrence", taxon12_last);
        StratigraphicRange sr13 = new StratigraphicRange();
        Taxon taxon13_first = new Taxon("t32_49");
        Taxon taxon13_last = new Taxon("t32_51");
        sr13.setInputValue("firstOccurrence", taxon13_first);
        sr13.setInputValue("lastOccurrence", taxon13_last);

        ArrayList<StratigraphicRange> sranges = new ArrayList<>();
        sranges.add(sr1);
        sranges.add(sr2);
        sranges.add(sr3);
        sranges.add(sr4);
        sranges.add(sr5);
        sranges.add(sr6);
        sranges.add(sr7);
        sranges.add(sr8);
        sranges.add(sr9);
        sranges.add(sr10);
        sranges.add(sr11);
        sranges.add(sr12);
        sranges.add(sr13);
        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.assignFrom(tree_initial);

        SRangesBirthDeathSkylineModel model = new SRangesBirthDeathSkylineModel();
        model.setInputValue("tree", tree);
        model.setInputValue("intervalTimes", new RealParameter("3 2 1.5 1 0"));
        model.setInputValue("origin", new RealParameter("4.0"));
        model.setInputValue("removalProbability", new RealParameter("0.0 0.0 0.0 0.0 0.0"));
        model.setInputValue("rho", new RealParameter("0.5"));
        model.setInputValue("netDiversification", new RealParameter("1.0 1.0 1.0 1.0 1.0"));
        model.setInputValue("turnOver", new RealParameter("0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333"));
        model.setInputValue("samplingProportion", new RealParameter("0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667"));

//        model.setInputValue("intervalTimes", new RealParameter("0"));
//        model.setInputValue("origin", new RealParameter("4.0"));
//        model.setInputValue("removalProbability", new RealParameter("0.0"));
//        model.setInputValue("rho", new RealParameter("0.5"));
//        model.setInputValue("netDiversification", new RealParameter("1.0"));
//        model.setInputValue("turnOver", new RealParameter("0.3333333333"));
//        model.setInputValue("samplingProportion", new RealParameter("0.1666666667"));
//        model.setInputValue("contemp", true);

        model.setInputValue("conditionOnSurvival", false);
        model.setInputValue("conditionOnRhoSampling", false);
        model.setInputValue("conditionOnRoot", false);

        model.initAndValidate();
        double ans = model.calculateLogP();
        System.out.println(ans);
        assertEquals(-283.60107898017947,ans , 1e-7);
    }
}
