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
        System.out.println("RHO");
        for (int i = 0; i < rho.length; i++) {
            System.out.print(rho[i] + " ");
        }

    }



//    public void birthExceedsDeath() {
//        Arrays.fill(birthExceedsDeath, false);
//    }


    public double Ai(int index){
//        System.out.println("Ai = "  + Ai(birth[index], death[index], psi[index]));
        return Ai(birth[index], death[index], psi[index]);
    }


    //KT interim: Correct in skyline
    //KT interim: Are the time indices off by one?
    public double q(double t) {
        int index = index(t);
        return q(t, index);
    }

    public double q(double t, int index) {
        return q(t, t_j(index), Ai(index), Bi(index));
    }

    public double q(double t, double t_i, double Ai, double Bi) {
        double v = Math.exp(-Ai * (t-t_i));
        return 4 * v / Math.pow((v*(1-Bi) + (1+Bi)), 2.0);
    }

    //KT interim: Are the time indices off by one?
    public double log_q(double t) {
        int index = index(t);
        return log_q(t, index);
    }

    public double log_q(double t, int index) {
        return log_q(index, t_j(index), t);

    }

    public double log_q(int index, double ti, double t) {
        // replacing Math.log( g(...) ) for better numerical stability
        return Math.log(4) - Ai(index) * (t - ti) - 2 * Math.log(Math.exp(-Ai(index) * (t - ti)) * (1 - Bi(index)) + (1 + Bi(index)));
    }

    public double q_tilde(double t){
        int index = index(t);
//        System.out.println("q_tilde at " + t);
//        System.out.println("t_j = " + t_j(t));
        return q_tilde(t, index);
    }

    public double q_tilde(double t, int index){
        return Math.sqrt(q(t,index)*Math.exp(-(birth[index]+death[index]+psi[index])*(t - t_j(index))));
    }

    public double log_q_tilde(double t){
        int index = index(t);
        return log_q_tilde(t, index);
    }

    public double log_q_tilde(double t, int index){
        return 0.5*(log_q(t,index)-(birth[index]+death[index]+psi[index])*(t-t_j(index)));
    }

    public double log_Q_tilde(double t){
    double result = 0;
        for (int k = 0; k < index(t); k++) {
        if (times[k] != t){
            result += log_q_tilde(times[k]);
        }
    }
        return result + log_q_tilde(t);
    }

    public double log_Q(double t){
        double result = 0;
        for (int k = 0; k < index(t); k++) {
            if (times[k] != t){
                result += log_q(times[k]);
            }
        }
        return result + log_q(t);
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
//                System.out.println("adjusting p where t = " + t + " and index = " + index);
                return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t, 0);
            }
        }
        return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t);
    }

    //KT interim: Set rho_j to zero in skyline model equation.
    public double Bi(int index){
//        System.out.println("calling single Bi at index " + index);
//        if (index == 0) {
//            System.out.println("Bi = " + Bi(birth[index], death[index], psi[index], 0, Ai(index), p_minus_1(index)));
//        }
        return Bi(birth[index], death[index], psi[index], 0, Ai(index), p_minus_1(index));
    }

    public double p_minus_1(int index) {
        if (index == 0){
//            System.out.println("p_minus_1 = " + (1-rho[totalIntervals-1]));
            return 1 - rho[totalIntervals-1];
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
            return 1 - rho[totalIntervals-1];
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
        System.out.println("node count: " + nodeCount);
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
            logP += log_Q(x0);
            System.out.println("conditioning 1 " + logP);
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return Double.NEGATIVE_INFINITY;
            } else {
                logP += log_Q(x1);
                System.out.println("logP " + logP);
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
        // integrate over fossils in the range. This seems to suggest that we take out the psi in the previous equations
        for (StratigraphicRange range : ((SRTree) tree).getSRanges()) {
            Node first = tree.getNode(range.getNodeNrs().get(0));
            List<Integer> fossils = range.getNodeNrs();
            int first_node_index = index(first.getHeight());
            double tFirst = first.getHeight();


            if (!range.isSingleFossilRange()) {
                Node lastNode = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size() - 1));
                double tLast = lastNode.getHeight();
                int last_node_index = index(tLast);

                System.out.println("SR branch from " + tFirst + " to " + tLast);
                logP += log_Q_tilde(tFirst) - log_Q(tFirst) ;
                logP += log_Q(tLast) - log_Q_tilde(tLast);
                System.out.println("minus q tilde " + tLast + " add q " + tLast);
                System.out.println("minus q " + tFirst + " add q tilde " + tFirst);



                if (first_node_index == last_node_index) {
                    logP += psi[first_node_index] * (tFirst - tLast); // TERM 3
                } else {
                    logP += psi[first_node_index] * (tFirst - times[first_node_index-1]); // TERM 13
                    logP += psi[last_node_index] * (times[last_node_index] - tLast); // TERM 14
                    for (int j =  last_node_index + 1; j < first_node_index; j++) {
                        logP += psi[j] * (times[j] - times[j-1]); // TERM 15
                    }
                }
            }

            Node ancestralLast = findAncestralRangeLastNode(first);
            // SR i is in I, i and a(i) are in the same vertical line in the graphical representation.
            if (ancestralLast != null) {
                double tAncestor = ancestralLast.getHeight();
                double tChild = first.getHeight();
                int ancestorIndex = index(tAncestor);
                int childIndex = index(tChild);
                double qsum = 1.0;
                for (int z = childIndex; z < ancestorIndex; z++) {
                    qsum *= q_tilde(times[z])/q(times[z]);
                }
                qsum = qsum*q_tilde(tAncestor)*q(tChild)/q_tilde(tChild)/q(tAncestor);
                qsum = 1 - qsum;
                logP += Math.log(qsum);

            }
        }

        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);
            double height = node.getHeight();
            Node parent = node.getParent();
            double parent_height = 0;
            int node_index = index(height);
                if (node.isLeaf()) {
                    if (!node.isDirectAncestor()) {
                        Node fossilParent = node.getParent();
                        if (node.getHeight() > 0.000000000005 || rho[totalIntervals-node_index-1] == 0.) {
                            logP += Math.log(p(height)); // - log_q(height); // TERM 4
                            if (((SRTree)tree).belongToSameSRange(i, fossilParent.getNr())) {
                                logP -= log_Q(height);

                            }
                            else {
                                logP -= log_Q(height);

                            }
                            logP += Math.log(psi[node_index]); // TERM 5
                        } else {
                            logP += Math.log(rho[totalIntervals-node_index-1]); //TERM 1
                        }
                    }

                } else {
                    if (node.isFake()) {
                        logP += Math.log(psi[node_index]);
                    } else {
                        logP += Math.log(birth[node_index]) + log_Q(height);
                    }
                }
            }
            System.out.println("END");
        return logP;
    }

}
