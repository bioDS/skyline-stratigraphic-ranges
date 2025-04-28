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
 * Edited by Kate Truman 2025
 */
public class SRangesBirthDeathSkylineModelTest extends TestCase {

    @Test
    public void testLikelihood() throws Exception {

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
        model.setInputValue("origin", new RealParameter("7.0"));
        model.setInputValue("removalProbability", new RealParameter("0.0 0.0 0.0 0.0 0.0 0.0 0.0"));
        model.setInputValue("rho", new RealParameter("0.5"));

        model.setInputValue("netDiversification", new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"));
        model.setInputValue("turnOver", new RealParameter("0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333 0.3333333333"));
        model.setInputValue("samplingProportion", new RealParameter("0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667 0.1666666667"));
        model.setInputValue("contemp", true);
        model.setInputValue("conditionOnSurvival", false);
        model.initAndValidate();
        double ans = model.calculateLogP();
        System.out.println(ans);
        assertEquals(-33.74668640318646, ans, 1e-9);


    }

    @Test
    public void testLikelihoodSingletonSR() throws Exception {

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
        model.setInputValue("rho", new RealParameter("0.7"));
        model.setInputValue("netDiversification", new RealParameter("0.25 0.5 0.09 0.5 0.5"));
        model.setInputValue("turnOver", new RealParameter("0.1666666667 0.2857142857 0.1 0.375 0.375"));
        model.setInputValue("samplingProportion", new RealParameter("0.6666666667 0.6 0.9756097561 0.7692307692 0.7692307692"));
        model.setInputValue("contemp", true);
        model.setInputValue("conditionOnSurvival", false);
//        model.setInputValue("reverseTimeArrays", true);


        model.initAndValidate();
        double ans = model.calculateLogP();
        System.out.println(ans);
        assertEquals(-24.68765411973161,ans , 1e-7);
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
        model.setInputValue("conditionOnSurvival", false);

        model.initAndValidate();
        double ans = model.calculateLogP();
        System.out.println(ans);
        assertEquals(-283.60107898017947,ans , 1e-7);
    }
}
