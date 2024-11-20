package sr.speciation;

import bdsky.evolution.speciation.BirthDeathSkylineModel;
import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import sa.evolution.speciation.SABirthDeathModel;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRTree;

import java.sql.Array;
import java.util.ArrayList;
import java.util.Arrays;
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

    }


    
    public void birthExceedsDeath() {
        Arrays.fill(birthExceedsDeath, false);
    }

    public double Ai(int index){
        return Ai(birth[index], death[index], psi[index]);
    }

    //KT interim: Set rho_j to zero in skyline model equation.
    public double Bi(int index){
        return Bi(birth[index], death[index], psi[index], rho[index], Ai(index), p[index]);
//        return (1 - 2(1 - p[index]))/Ai(index);
    }


    //KT interim: Correct in skyline
    //KT interim: Are the time indices off by one?
    public double q(double t) {
        int index = index(t);
        return q(t, Ai(index), Bi(index));
    }

    public double q(double t, double Ai, double Bi) {
        double v = Math.exp(-Ai * t);
        return 4 * v / Math.pow(v*(1-Bi) + (1+Bi), 2.0);
    }

    //KT interim: Are the time indices off by one?
    public double log_q(double t) {
        int index = index(t);
        return log_q(index, times[index], t);
    }

    public double q_tilde(double t){
        int index = index(t);
        return Math.sqrt(q[index]*Math.exp(-(birth[index]+death[index]+psi[index]))*(t - times[index]));
    }

    public double log_q_tilde(double t){
        int index = index(t);
        return 0.5*(log_q(t)-(birth[index]+death[index]+psi[index])*(t-times[index]));
    }

    public double p(double t) {
        int index = index(t);
        return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t);
    }

    //KT: equation from FBD skyline document
    public double p(double birth, double death, double psi, double A, double B, double t){
        int index = index(t);
        return (birth + death + psi - A*((1+ B)-(1-B)*Math.exp(-A*(t - times[index])))/((1+ B)-(1-B)*Math.exp(-A*(t - times[index])))/2*birth);
    }

    public double log_p0s(double t) {
        int index = index(t);
        return log_p0s(t, Ai(index), Bi(index), index);
    }

    public double log_p0s(double t, double Ai, double Bi, int index) {
        double p0 = (birth[index] + death[index] + psi[index] - Ai * ((1 + Bi) - Math.exp(-Ai * t) * (1 - Bi)) / ((1 + Bi) + Math.exp(-Bi * t) * (1 - Bi))) / (2 * birth[index]);
        return Math.log(r[index] + (1 - r[index]) * p0);
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
        // for node in nodes:
        // if node isDirectAncestor (SA):


        SRTree tree = (SRTree) treeInput.get();
        int nodeCount = tree.getNodeCount();
        //TODO: is a replacement update needed?
//        updateParameters();
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
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return Double.NEGATIVE_INFINITY;
            } else {
                logP += log_q(x1);
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
        //check for multi-time branch for later.
        //if(findIndex(parent_height) == node_index){}

        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);
            double height = node.getHeight();
            Node parent = node.getParent();
            double parent_height = -1;
            int node_index = index(height);
            int sub_parent = 0;
            double skyline_parent = -1;
            if (node.getChildCount() == 2){
                logP += birth[node_index];
            }
            else {
                if(!node.isFake()){
                    parent_height = parent.getHeight();
                    skyline_parent = parent_height;
                    if (times[node_index-1] < parent_height && times[node_index-1] > height){
                        skyline_parent = times[node_index-1];
                    }
                    // Not SR branch
                    // TODO: are q_tilde and q the right way around for time intervals?
                    if(tree.getRangeOfNode(node) == null){
                        logP += log_q(parent_height)-log_q(height);
                        if (skyline_parent != parent_height){
                            logP += Math.log(q(skyline_parent));
                        }
                    }
                    else{
                        logP += log_q_tilde(parent_height)-log_q_tilde(height);
                        if (skyline_parent != parent_height){
                            logP += Math.log(q_tilde(skyline_parent));
                        }
                    }
                }

                if (node.isLeaf()) {
                    if  (!node.isDirectAncestor())  {
                        Node fossilParent = node.getParent();
                        if (node.getHeight() > 0.000000000005 || rho[node_index] == 0.) {
                            logP += Math.log(p(height));
    //                            //SR y_i

    //                        if (((SRTree)tree).belongToSameSRange(i, fossilParent.getNr())) {
    //                            logP += Math.log(psi[node_index]) - log_q_tilde(node.getHeight()) + log_p0s(node.getHeight());
    //                        } else {
    //                            // SR o_i
    //                            logP += Math.log(psi[node_index]) - log_q(node.getHeight()) + log_p0s(node.getHeight());
    //                        }
                        } else {
                            logP += Math.log(rho[node_index]); //chi_i ln(rho) term
                        }
                    }
                } else {
                    if (node.isFake()) {
    //                    logP += Math.log(psi[node_index]);
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
                    }
                }
            }

            // integrate over fossils in the range. This seems to suggest that we take out the psi in the previous equations
            for (StratigraphicRange range:((SRTree)tree).getSRanges()) {
                Node first =  tree.getNode(range.getNodeNrs().get(0));
                List<Integer> fossils = range.getNodeNrs()
                int first_node_index = index(first.getHeight());
                logP += Math.log(psi[index(first.getHeight())]);
                if (fossils.size() != 1) {
                    logP += Math.log(psi[index(tree.getNode(fossils.get(fossils.size()) - 1).getHeight())]);
                }

            }
                if (!range.isSingleFossilRange()) {
                    double tFirst =first.getHeight();
                    double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size()-1)).getHeight();
                    int last_node_index = index(tLast);
                    if (first_node_index == last_node_index) {
                        logP += psi[first_node_index] * (tFirst - tLast);
                    }
                    else {
                        logP += psi[first_node_index] * (tFirst - times[first_node_index]);
                        logP += psi[last_node_index] * (times[last_node_index] -tLast);
                        for (int j = first_node_index + 1; j < last_node_index; j++) {
                            logP += psi[j]*(times[j-1] -times[j]);
    //                        logP += Math.log(q(times[j-1]));
                        }
                    }
                }
                Node ancestralLast = findAncestralRangeLastNode(first);
                // SR i is in I, i and a(i) are in the same vertical line in the graphical representation.
                 if (ancestralLast != null) {
                     double tOld = ancestralLast.getHeight();
                     double tYoung = first.getHeight();
                     int oldIndex = index(tOld);
                     int youngIndex = index(tYoung);
                     double qsum;
                     if (oldIndex == youngIndex) {
                         logP += Math.log(1 - q(tYoung) / q_tilde(tYoung) * q_tilde(tOld) / q(tOld));
                     } else {
                         qsum = 1 - q(times[oldIndex]) / q(tOld) * q_tilde(tOld) / q_tilde(times[oldIndex]);
                         qsum += 1 - q(youngIndex) / q(times[youngIndex]) * q_tilde(times[youngIndex - 1]) * q_tilde(times[youngIndex - 1]) / q_tilde(youngIndex);
                         for (int j = oldIndex + 1; j < youngIndex; j++) {
                             qsum += 1 - q(times[j]) / q(times[j - 1]) * q_tilde(times[j - 1] / q_tilde(times[j]));
                         }
                         logP += Math.log(qsum);
                     }
                 }
             }
        return logP;
    }

}
