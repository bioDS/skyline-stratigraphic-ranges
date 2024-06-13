package sr.speciation;

import bdsky.evolution.speciation.BirthDeathSkylineModel;
import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import sa.evolution.speciation.SABirthDeathModel;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRSkylineTree;

import java.util.ArrayList;
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
public class SRangesBirthDeathSkylineModel extends SABirthDeathModel {

    // Replace constant values from parent class.

    //TODO: Why does BirthDeathSkylineModel class only use singular values but comments describe multiple?
    public Input<List<RealParameter>> birthRateInput = new Input<List<RealParameter>>("birthRates","Birth rates",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);
    public Input<List<RealParameter>> deathRateInput = new Input<List<RealParameter>>("deathRates","Death rates",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);
    public Input<List<RealParameter>> samplingRateInput = new Input<List<RealParameter>>("samplingRates","Sampling rates",
            new ArrayList<RealParameter>(), Input.Validate.REQUIRED);

    public Input<List<RealParameter>> diversificationRateInput = new Input<List<RealParameter>>("diversificationRates",
            "Net diversification rates. Birth rate - death rate", Input.Validate.XOR, birthRateInput);
    public Input<List<RealParameter>> turnoverInput = new Input<List<RealParameter>>("turnovers",
            "Turnovers. Death rate / birth rate", Input.Validate.XOR, deathRateInput);
    public Input<List<RealParameter>> samplingProportionInput = new Input<List<RealParameter>>("samplingProportions",
            "The probability of sampling prior to death. Sampling rate /(sampling rate + death rate)", Input.Validate.XOR, samplingRateInput);

    // Replace constant values from parent class.
    protected Double[] lambda;
    protected Double[] mu;
    protected Double[] psi;
    protected Double[] p0;


    protected BirthDeathSkylineModel skyline = new BirthDeathSkylineModel();

    public double Ai(int index){
        return skyline.Ai(lambda[index], mu[index], psi[index]);
    }

    public double Bi(int index){
        return skyline.Bi(lambda[index], mu[index], psi[index], rho, Ai(index), p0[index]);
//        return (1 - 2(1 - p[index]))/Ai(index);
    }

    @Override
    public double log_q(int index, double ti, double t) {
        return skyline.log_q(index, ti, t);
    }

    public double q_tilde(int index, double ti, double t){
        return Math.sqrt(q[index]*Math.exp(-(birth[index]+death[index]+psi[index]))*(t - ti));
    }

    // p0 fine as is in BirthDeathSkylineModel?

//    private double log_q_tilde(double t, double c1, double c2) {
//        return 0.5*(-t*(lambda + mu + psi) + log_q(t,c1,c2));
//    }

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
        SRTree tree = (SRTree) treeInput.get();
        int nodeCount = tree.getNodeCount();
        updateParameters();
        if (lambdaExceedsMu && lambda <= mu) {
            return Double.NEGATIVE_INFINITY;
        }

        if (lambda < 0 || mu < 0 || psi < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        double x0 = origin;
        double x1=tree.getRoot().getHeight();

        if (x0 < x1 ) {
            return Double.NEGATIVE_INFINITY;
        }

        if (!conditionOnRootInput.get()){
            logP = log_q(x0, c1, c2);
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return Double.NEGATIVE_INFINITY;
            } else {
                logP = log_q(x1, c1, c2);
            }
        }

        if (conditionOnSamplingInput.get()) {
            logP -= log_oneMinusP0(x0, c1, c2);
        }

        if (conditionOnRhoSamplingInput.get()) {
            if (conditionOnRootInput.get()) {
                logP -= Math.log(lambda) + log_oneMinusP0Hat(x1, c1, c2)+ log_oneMinusP0Hat(x1, c1, c2);
            }  else {
                logP -= log_oneMinusP0Hat(x0, c1, c2);
            }
        }

        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);
            if (node.isLeaf()) {
                if  (!node.isDirectAncestor())  {
                    Node fossilParent = node.getParent();
                    if (node.getHeight() > 0.000000000005 || rho == 0.) {

                        if (((SRTree)tree).belongToSameSRange(i, fossilParent.getNr())) {
                            logP += Math.log(psi) - log_q_tilde(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                        } else {
                            logP += Math.log(psi) - log_q(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                        }
                    } else {
                        logP += Math.log(rho);
                    }
                }
            } else {
                if (node.isFake()) {
                    logP += Math.log(psi);
                    Node parent = node.getParent();
                    Node child = node.getNonDirectAncestorChild();
                    Node DAchild = node.getDirectAncestorChild();
                    if (parent != null && ((SRTree)tree).belongToSameSRange(parent.getNr(),DAchild.getNr())) {
                        logP += - log_q_tilde(node.getHeight(), c1, c2) + log_q(node.getHeight(), c1, c2);
                    }
                    if (child != null && ((SRTree)tree).belongToSameSRange(i,child.getNr())) {
                        logP += - log_q(node.getHeight(), c1, c2) +  log_q_tilde(node.getHeight(), c1, c2);
                    }
                } else {
                    logP += Math.log(lambda) + log_q(node.getHeight(), c1, c2);
                }
            }
        }

        // integrate over fossils in the range. This seems to suggest that we take out the psi in the previous equations
        for (StratigraphicRange range:((SRTree)tree).getSRanges()) {
            Node first =  tree.getNode(range.getNodeNrs().get(0));
            if (!range.isSingleFossilRange()) {
                double tFirst =first.getHeight();
                double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size()-1)).getHeight();
                logP += psi*(tFirst - tLast);
            }
            Node ancestralLast = findAncestralRangeLastNode(first);
            if (ancestralLast != null) {
                double tOld = ancestralLast.getHeight();
                double tYoung = first.getHeight();
                logP += Math.log(1-q(tYoung, c1, c2)/q_tilde(tYoung, c1, c2)*q_tilde(tOld, c1, c2)/q(tOld, c1, c2));
            }
        }
        return logP;
    }

}
