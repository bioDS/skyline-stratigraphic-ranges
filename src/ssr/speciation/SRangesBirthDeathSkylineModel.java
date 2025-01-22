package ssr.speciation;

import bdsky.evolution.speciation.BirthDeathSkylineModel;
import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRTree;

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

    protected Double[] p;
    protected Double[] q;
    protected int l;
    protected Boolean birthExceedsDeath[] = new Boolean[l];

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
        super.transformParameters_d_r_s();

    }

    public double Ai(int index){
        return Ai(birth[index], death[index], psi[index]);
    }

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

    public double log_q(double t) {
        int index = index(t);
        return log_q(t, index);
    }

    public double log_q(double t, int index) {
        return log_q(index, t_j(index), t);

    }

    public double log_q(int index, double ti, double t) {
        return Math.log(4) - Ai(index) * (t - ti) - 2 * Math.log(Math.exp(-Ai(index) * (t - ti)) * (1 - Bi(index)) + (1 + Bi(index)));
    }

    public double q_tilde(double t){
        int index = index(t);
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

    public double log_q_comb_tilde(double t){
    double result = 0;
        for (int k = 0; k < index(t); k++) {
        if (times[k] != t){
            result += log_q_tilde(times[k]);
        }
    }
        return result + log_q_tilde(t);
    }

    public double log_q_comb(double t){
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
                return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t, t_j(index));
            }
            else {
                return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t, 0);
            }
        }
        return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t);
    }

    public double Bi(int index){
        return Bi(birth[index], death[index], psi[index], 0, Ai(index), p_minus_1(index));
    }

    public double p_minus_1(int index) {
        if (index == 0){
            return 1 - rho[totalIntervals-1];
        }
        index = index-1;
        double t = times[index];
        return p(birth[index], death[index], psi[index], Ai(index), Bi(index), t, t_j(index));
    }

    public double p(double birth, double death, double psi, double A, double B, double t){
        return p(birth, death, psi, A, B,  t, t_j(t));
    }

    public double p(double birth, double death, double psi, double A, double B, double t, double t_j) {
        if (t == 0){
            return 1 - rho[totalIntervals-1];
        }
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
            logP += log_q_comb(x0);
            System.out.println("conditioning 1 " + logP);
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return Double.NEGATIVE_INFINITY;
            } else {
                logP += log_q_comb(x1);
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
        for (StratigraphicRange range : tree.getSRanges()) {
            Node first = tree.getNode(range.getNodeNrs().get(0));
            double tFirst = first.getHeight();
            int first_node_index = index(tFirst);

            if (!range.isSingleFossilRange()) {
                Node lastNode = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size() - 1));
                double tLast = lastNode.getHeight();
                int last_node_index = index(tLast);

                logP += log_q_comb_tilde(tFirst) - log_q_comb(tFirst);
                logP += log_q_comb(tLast) - log_q_comb_tilde(tLast);

                if (first_node_index == last_node_index) {
                    logP += psi[first_node_index] * (tFirst - tLast);
                } else {
                    logP += psi[first_node_index] * (tFirst - times[first_node_index-1]);
                    logP += psi[last_node_index] * (times[last_node_index] - tLast);
                    for (int j =  last_node_index + 1; j < first_node_index; j++) {
                        logP += psi[j] * (times[j] - times[j-1]);
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
                logP += Math.log(1-qsum);
            }
        }

        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);
            double height = node.getHeight();
            int node_index = index(height);
                if (node.isLeaf()) {
                    if (!node.isDirectAncestor()) {
                        Node fossilParent = node.getParent();
                        if (height > 0.000000000005 || rho[totalIntervals-node_index-1] == 0.) {
                            logP += Math.log(p(height));
                            if (tree.belongToSameSRange(i, fossilParent.getNr())) {
                                logP -= log_q_comb(height);

                            }
                            else {
                                logP -= log_q_comb(height);

                            }
                            logP += Math.log(psi[node_index]);
                        } else {
                            logP += Math.log(rho[totalIntervals-node_index-1]);
                        }
                    }

                } else {
                    if (node.isFake()) {
                        logP += Math.log(psi[node_index]);
                    } else {
                        logP += Math.log(birth[node_index]) + log_q_comb(height);
                    }
                }
            }
        return logP;
    }

}
