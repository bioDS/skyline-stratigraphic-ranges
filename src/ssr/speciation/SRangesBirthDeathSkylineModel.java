package ssr.speciation;

import bdsky.evolution.speciation.BirthDeathSkylineModel;
import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRTree;

import java.util.List;

/**
 * @author Alexandra Gavryushkina
 * @author Ugne Stolz
 * @author Kate Truman
 */


@Description("A variant of the fossilized birth-death model under budding (asymmetric) speciation with stratigraphic ranges")
@Citation("Stadler T, Gavryushkina A, Warnock RCM, Drummond AJ, Heath TA (2017) \n" +
        "The fossilized birth-death model under different types of speciation")
@Citation("Gavryushkina A, Warnock RCM, Drummond AJ, Heath TA, Stadler T (2017) \n" +
        "Bayesian total-evidence dating under the fossilized birth-death model with stratigraphic ranges.")
public class SRangesBirthDeathSkylineModel extends BirthDeathSkylineModel {

//    // the interval times
//    public Input<RealParameter> timesInput =
//            new Input<RealParameter>("times", "The times t_i specifying when diversification rate changes occur", (RealParameter) null);

    

    // Replace constant values from parent class.

    //KT interim: Why does BirthDeathSkylineModel class only use singular values but comments describe multiple?
//    public Input<List<RealParameter>> birthRateInput = new Input<List<RealParameter>>("birthRates","Birth rates",
//            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);
//    public Input<List<RealParameter>> deathRateInput = new Input<List<RealParameter>>("deathRates","Death rates",
//            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);
//    public Input<List<RealParameter>> samplingRateInput = new Input<List<RealParameter>>("samplingRates","Sampling rates",
//            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);
//
//    public Input<List<RealParameter>> diversificationRateInput = new Input<List<RealParameter>>("diversificationRates",
//            "Net diversification rates. Birth rate - death rate", Input.Validate.XOR, birthRateInput);
//    public Input<List<RealParameter>> turnoverInput = new Input<List<RealParameter>>("turnovers",
//            "Turnovers. Death rate / birth rate", Input.Validate.XOR, deathRateInput);
//    public Input<List<RealParameter>> samplingProportionInput = new Input<List<RealParameter>>("samplingProportions",
//            "The probability of sampling prior to death. Sampling rate /(sampling rate + death rate)", Input.Validate.XOR, samplingRateInput);

    protected Double[] p;
    protected Double[] q;
    // Number of skyline intervals
    protected int l;
    protected Boolean birthExceedsDeath[] = new Boolean[l];
//    protected List<Double> times = new ArrayList<Double>();
    Boolean timesRelative = false;




    // protected List<Boolean> birthExceedsDeath = new ArrayList<Boolean>();

//    @Override
//    public void setInputValue(String name, Object value){
////        new List<String> skylineParams = new ArrayList<>(Arrays.asList("",""))
//            try {
//                skyline.setInputValue(name, value);
//            }
//            catch (RuntimeException e){
//                super.setInputValue(name, value);
//        }
//    }

    //TODO: what needs adding to this?
    public void initAndValidate() {
        int timeInputs = 0;
        if (intervalTimes.get() != null) {
            timeInputs += 1;
        }
        if (birthRateChangeTimesInput.get() != null) {
            timeInputs += 1;
            intervalTimes = birthRateChangeTimesInput;
        }
        if (deathRateChangeTimesInput.get() != null) {
            timeInputs += 1;
            intervalTimes = deathRateChangeTimesInput;

        }
        if (samplingRateChangeTimesInput.get() != null) {
            timeInputs += 1;
            intervalTimes = samplingRateChangeTimesInput;
        }
        if (removalProbabilityChangeTimesInput.get() != null) {
            timeInputs += 1;
            intervalTimes = removalProbabilityChangeTimesInput;
        }

        if (timeInputs > 1){
            throw new RuntimeException("More than one set of diversification rate change times specified." +
                    "Please specify only the overall time intervals, or either the birth, death, sampling or removal probability change times" +
                    "as a proxy for this.");
        }
        super.initAndValidate();
//        super.updateRatesAndTimes();
        System.out.println("did super init");
        super.transformParameters_d_r_s();
//        for (int i = 0; i < birth.length; i++) {
//            System.out.print(birth[i] + " ");
//        }
//        System.out.println();
//        for (int i = 0; i < death.length; i++) {
//            System.out.print(death[i] + " ");
//        }
//        System.out.println();
//        for (int i = 0; i < psi.length; i++) {
//            System.out.print(psi[i] + " ");
//        }

    }


    
//    public void birthExceedsDeath() {
//        Arrays.fill(birthExceedsDeath, false);
//    }


    public double Ai(int index){
//        System.out.println("calling single Ai");
        System.out.println("Ai = "  + Ai(birth[index], death[index], psi[index]));
        return Ai(birth[index], death[index], psi[index]);
    }




    //KT interim: Correct in skyline
    //KT interim: Are the time indices off by one?
    public double q(double t) {
        int index = index(t);
        return q(t, Ai(index), Bi(index));
    }

    public double q(double t, double Ai, double Bi) {
        int index = index(t);
        double v = Math.exp(-Ai * (t-t_j(t)));
        System.out.println("q at " + t + " is " + (4 * v / Math.pow((v*(1-Bi) + (1+Bi)), 2.0)));
        System.out.println("q values t " + t + " " + t_j(t));
        return 4 * v / Math.pow((v*(1-Bi) + (1+Bi)), 2.0);
    }

    //KT interim: Are the time indices off by one?
    public double log_q(double t) {
        int index = index(t);
        return log_q(index, t_j(t), t);
    }

    public double q_tilde(double t){
        int index = index(t);
        System.out.println("q_tilde at " + t);
        System.out.println("t_j = " + t_j(t));
        return Math.sqrt(q(t)*Math.exp(-(birth[index]+death[index]+psi[index])*(t - t_j(t))));
    }

    public double log_q_tilde(double t){
        int index = index(t);
        return 0.5*(log_q(t)-(birth[index]+death[index]+psi[index])*(t-t_j(t)));
    }

    public double t_j(double t){
        int index = index(t);
        if(t == times[index]){
            return t;
        }
        if (index != 0){
            return times[index - 1];
        }
        return 0;
    }

    public double t_j(int index){
    if (index == 0){
        return 0;
    }
    return times[index-1];
    }

    public double p(double t) {
    int index = index(t);
        if (times[index] == t){
            if (index != 0){
//                index = index - 1;

                return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t, t_j(index));
            }
            else {
                System.out.println("adjusting p where t = " + t + " and index = " + index);
                return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t, 0);
            }
        }
        return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t);
    }

    //KT interim: Set rho_j to zero in skyline model equation.
    public double Bi(int index){
//        System.out.println("calling single Bi at index " + index);
        if (index == 0) {
            System.out.println("Bi = " + Bi(birth[index], death[index], psi[index], rho[index], Ai(index), p_minus_1(index)));
        }
        return Bi(birth[index], death[index], psi[index], rho[index], Ai(index), p_minus_1(index));
    }

    public double p_minus_1(int index) {
        if (index == 0){
            System.out.println("p_minus_1 = " + (1-rho[0]));
            return 1 - rho[0];
        }
//        System.out.println("p_add_1 at index " + index);
        index = index-1;
        double t = times[index];
//        System.out.println("t = " + t + " t_j " + t_j(index));
//        System.out.println("index" + index);
        return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t, t_j(index));
    }



    //KT: equation from FBD skyline document
    public double p(double birth, double death, double psi, double A, double B, double t){
        return p(birth, death, psi, A, B,  t, t_j(t));
    }

    public double p(double birth, double death, double psi, double A, double B, double t, double t_j) {
        if (t == 0){
            return 1 - rho[0];
        }
//        System.out.println("t_j = " + t_j(t));
        double mid = A * (1 + B - (1 - B) * Math.exp(-A * (t - t_j))) / (1 + B + (1 - B) * Math.exp(-A * (t - t_j)));
        return (birth + death + psi - mid)/(2*birth);

    }



        private Node findAncestralRangeLastNode(Node node) {
        Node parent = node.getParent();
        if (node.isDirectAncestor()){
            parent = parent.getParent();
            node = node.getParent();
        }
        if (parent == null) {
            return parent;
        } else {
            if (parent.isFake()) {
                return parent;
            } else if (parent.getChild(0) == node) {
                return findAncestralRangeLastNode(parent);
            } else {
                return null;
            }
        }
    }



    @Override
    public double calculateLogP()
    {
        double logP = 0;
        SRTree tree = (SRTree) treeInput.get();
        int nodeCount = tree.getNodeCount();
        preCalculation(tree);
        for(int i = 0; i < l; i++) {
            if (birthExceedsDeath[i] && birth[i] <= death[i]) {
                return Double.NEGATIVE_INFINITY;
            }
            if (birth[i] < 0 || death[i] < 0 || psi[i] < 0) {
                return Double.NEGATIVE_INFINITY;
            }
        }



        double x0 = origin.get().getArrayValue();
        double x1=tree.getRoot().getHeight();

        if (x0 < x1 ) {
            return Double.NEGATIVE_INFINITY;
        }

        if (!conditionOnRootInput.get()){
            logP += log_q(x0);
            System.out.println("logP" + logP);
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return Double.NEGATIVE_INFINITY;
            } else {
                logP += log_q(x1);
                System.out.println("logP" + logP);
            }
        }

       //TODO: reinstate with corrections?
//        if (conditionOnSamplingInput.get()) {
//            logP -= log_oneMinusP0(x0, c1, c2);
//        }
//
//        if (conditionOnRhoSamplingInput.get()) {
//            if (conditionOnRootInput.get()) {
//                logP -= Math.log(lambda) + log_oneMinusP0Hat(x1, c1, c2)+ log_oneMinusP0Hat(x1, c1, c2);
//            }  else {
//                logP -= log_oneMinusP0Hat(x0, c1, c2);
//            }
//        }
        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);
            double height = node.getHeight();
            Node parent = node.getParent();
            double parent_height = -1;
            int node_index = index(height);
            double skyline_parent = -1;
            if (node.getChildCount() == 2) {
                logP += birth[node_index];  // TERM 2
                System.out.println("logP @ 2 =" + logP);

            } else {
                if (!node.isFake()) {
                    parent_height = parent.getHeight();
                    skyline_parent = parent_height;
                    if ((node_index < totalIntervals - 1) && (times[node_index + 1] < parent_height && times[node_index + 1] > height)) {
                        skyline_parent = times[node_index + 1];
                    }
                    // Not SR branch
                    // TODO: are q_tilde and q the right way around for time intervals?
                    if (tree.getRangeOfNode(node) == null) {
                        logP += log_q(parent_height) - log_q(height); // TERM 8
                        System.out.println("logP @ 8 =" + logP);

                        if (skyline_parent != parent_height) {
                            logP += Math.log(q(skyline_parent)); // TERM 6
                            System.out.println("logP @ 6a =" + logP);

                        }
                    } else {
                        logP += log_q_tilde(parent_height) - log_q_tilde(height);
                        System.out.println("logP" + logP);
                        System.out.println("logP @ 7 =" + logP);
                        if (skyline_parent != parent_height) {
                            logP += Math.log(q_tilde(skyline_parent));  // TERM 7
                            System.out.println("logP @ 6b =" + logP);
                        }
                    }
                }

                if (node.isLeaf()) {
                    if (!node.isDirectAncestor()) {
                        if (node.getHeight() > 0.000000000005 || rho[node_index] == 0.) {
                            logP += Math.log(p(height)); // TERM 4
                            System.out.println("adding p at " + height);
                            System.out.println("logP @ 4 =" + logP);
                            //                            //SR y_i

                            //                        if (((SRTree)tree).belongToSameSRange(i, fossilParent.getNr())) {
                            //                            logP += Math.log(psi[node_index]) - log_q_tilde(node.getHeight()) + log_p0s(node.getHeight());
                            //                        } else {
                            //                            // SR o_i
                            //                            logP += Math.log(psi[node_index]) - log_q(node.getHeight()) + log_p0s(node.getHeight());
                            //                        }
                        } else {
                            logP += Math.log(rho[node_index]); //TERM 1
                            System.out.println("logP @ 1 =" + logP);
                            System.out.println("adding rho");
                        }
                    }
                } else {
                    if (node.isFake()) {
                                            logP += Math.log(psi[node_index]);
                                            System.out.println("adding psi at " + height);
                        parent = node.getParent();
                        //                    Node child = node.getNonDirectAncestorChild();
                        //                    Node DAchild = node.getDirectAncestorChild();
                        //                    if (parent != null && ((SRTree)tree).belongToSameSRange(parent.getNr(),DAchild.getNr())) {
                        //                        logP += - log_q_tilde(node.getHeight()) + log_q(node.getHeight());
                        //                    }
                        //                    if (child != null && ((SRTree)tree).belongToSameSRange(i,child.getNr())) {
                        //                        logP += - log_q(node.getHeight()) +  log_q_tilde(node.getHeight());
                        //                    }
                    } else {
                        // speciation event
                        //                    logP += Math.log(birth[node_index]) + log_q(node.getHeight());
                        logP += Math.log(birth[node_index]);
                        System.out.println("logP @ unknown =" + logP);
                    }
                }
            }
        }

            // integrate over fossils in the range. This seems to suggest that we take out the psi in the previous equations
            for (StratigraphicRange range : ((SRTree) tree).getSRanges()) {
                Node first = tree.getNode(range.getNodeNrs().get(0));
                List<Integer> fossils = range.getNodeNrs();
                int first_node_index = index(first.getHeight());
//                logP += Math.log(psi[index(first.getHeight())]); // TERM 5
//                System.out.println("adding psi at " + first.getHeight());
                System.out.println("logP @ 5a =" + logP);

                if (fossils.size() != 1) {
//                    logP += Math.log(psi[index(tree.getNode(fossils.get(fossils.size()-1)).getHeight())]);
//                    System.out.println("adding psi at " + tree.getNode(fossils.get(fossils.size()-1)).getHeight());
                    System.out.println("logP @ 5b =" + logP);
                }


                if (!range.isSingleFossilRange()) {
                    double tFirst = first.getHeight();
                    double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size() - 1)).getHeight();
                    int last_node_index = index(tLast);
                    System.out.println("first node index: " + first_node_index + ", last node index: " + last_node_index);
                    System.out.println("Stratigraphic range from: " + tFirst + " to : " + tLast);
                    if (first_node_index == last_node_index) {
                        logP += psi[first_node_index] * (tFirst - tLast); // TERM 3
                        System.out.println("adding psi difference");
                        System.out.println("logP @ 3 =" + logP);
                    } else {
                        logP += psi[first_node_index] * (tFirst - times[first_node_index]); // TERM 13
                        logP += psi[last_node_index] * (times[last_node_index] - tLast); // TERM 14
                        System.out.println("adding psi difference");
                        System.out.println("logP @ 14 =" + logP);
                        for (int j =  last_node_index + 1; j < first_node_index; j++) {
                            logP += psi[j] * (times[j - 1] - times[j]); // TERM 15
                            System.out.println("logP @ 15 =" + logP);
                            //                        logP += Math.log(q(times[j-1]));
                        }
                    }
                }
                Node ancestralLast = findAncestralRangeLastNode(first);
                // SR i is in I, i and a(i) are in the same vertical line in the graphical representation.
                if (ancestralLast != null) {
                    System.out.println("ancestralLast is " + ancestralLast.getHeight() + " first = " + first.getHeight());

                    double tAncestor = ancestralLast.getHeight();
                    double tChild = first.getHeight();
                    int ancestorIndex = index(tAncestor);
                    int childIndex = index(tChild);
                    System.out.println("ancestor index: " + ancestorIndex);
                    System.out.println("child index: " + childIndex);
                    System.out.println("tAncestor: " + tAncestor);
                    System.out.println("tChild: " + tChild);
                    double qsum;
                    if (ancestorIndex == childIndex) {
                        logP += Math.log(1 - q(tChild) / q_tilde(tChild) * q_tilde(tAncestor) / q(tAncestor)); // TERM 9
                        System.out.println("logP @ 9 =" + logP);
                    } else {
                        System.out.println("tAncestor: " + tAncestor + "ancestorIndex: " + ancestorIndex + "time " + times[ancestorIndex]);
                       System.out.println("here:" + q(times[ancestorIndex-1]) + " " + q(tAncestor) + " "+ q_tilde(tAncestor) + " " + q_tilde(times[ancestorIndex-1]));
                       System.out.println("index of time 5 is " + index(5));
                        qsum = 1 - q(times[ancestorIndex-1]) / q(tAncestor) * q_tilde(tAncestor) / q_tilde(times[ancestorIndex-1]); // TERM 10
                        System.out.println("qsum a = " + qsum);
                        qsum += 1 - q(tChild) / q(times[childIndex]) * q_tilde(times[childIndex]) / q_tilde(tChild); // TERM 11
                        System.out.println("qsum b = " + (1 - q(tChild) / q(times[childIndex]) * q_tilde(times[childIndex]) / q_tilde(tChild)) + " " + childIndex + " " + tChild + " " + q_tilde(tChild));
                       for (int j = childIndex + 1; j < ancestorIndex; j++) {
                            qsum += 1 - q(times[j]) / q(times[j + 1]) * q_tilde(times[j + 1] / q_tilde(times[j])); // TERM 12
                            System.out.println("qsum c = " + (1 - q(times[j]) / q(times[j + 1]) * q_tilde(times[j + 1] / q_tilde(times[j]))));

                        }
                        System.out.println("qsum d = " + qsum);

                        logP += Math.log(qsum);
                        System.out.println("logP @ 12 =" + logP);
                    }
                }
            }
            System.out.println("END");
            System.out.println(p(2));
        return logP;
    }

}
