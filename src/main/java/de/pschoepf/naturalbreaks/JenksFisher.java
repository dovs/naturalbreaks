/*
 * Port of Jenks/Fisher breaks originally created in C by Maarten Hilferink.
 *
 * Copyright (C) {2015}  {Philipp Schoepf}
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/
package de.pschoepf.naturalbreaks;

import lombok.Data;
import lombok.RequiredArgsConstructor;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Basic port of original C code from Maarten Hilferink.
 * All credits this fantastic work go to him for.
 *
 * @author Philipp Schöpf
 */
public class JenksFisher {

    private List<ValueCountPair> cumulValues;
    private int numValues;
    private int numBreaks;
    private int bufferSize;
    private double[] previousSSM;
    private double[] currentSSM;
    private int[] classBreaks;
    private int classBreaksIndex;
    private int completedRows;

    /**
     * Constructor that initializes main variables used in fisher calculation of natural breaks.
     *
     * @param vcpc Ordered list of pairs of values to occurrence counts.
     * @param k    Number of breaks to find.
     */
    private JenksFisher(final int k, List<ValueCountPair> vcpc) {
        cumulValues = new ArrayList<>();
        numValues = vcpc.size();
        numBreaks = k;
        bufferSize = (vcpc.size() - (k - 1));
        previousSSM = new double[bufferSize];
        currentSSM = new double[bufferSize];
        classBreaks = new int[bufferSize * (numBreaks - 1)];
        classBreaksIndex = 0;
        completedRows = 0;
        double cwv = 0.0;
        long cw = 0;

        ValueCountPair currPair;

        for (int i = 0; i != numValues; ++i) {
            currPair = vcpc.get(i);
            assert (i == 0 || currPair.getValue() >= vcpc.get(i - 1).getValue()); // PRECONDITION: the value sequence must be strictly increasing
            long w = currPair.getCount();
            cw = Math.addExact(cw, w);
            cwv += w * currPair.getValue();
            cumulValues.add(new ValueCountPair(cwv, cw));
            if (i < bufferSize) {
                previousSSM[i] = cwv * cwv / cw; // prepare sum of squared means for first class. Last (k-1) values are omitted
            }
        }
    }

    /**
     * Gets sum of weighs for elements with index b..e.
     *
     * @param b index of begin element
     * @param e index of end element
     * @return sum of weights.
     */
    private long getSumOfWeights(int b, int e) {
        assert (b != 0);    // First element always belongs to class 0, thus queries should never include it.
        assert (b <= e);
        assert (e < numValues);

        long res = cumulValues.get(e).getCount();
        res -= cumulValues.get(b - 1).getCount();
        return res;
    }

    /**
     * Gets sum of weighed values for elements with index b..e
     *
     * @param b index of begin element
     * @param e index of end element
     * @return the cumul. sum of the values*weight
     */
    private double getSumOfWeightedValues(int b, int e) {
        assert (b != 0);
        assert (b <= e);
        assert (e < numValues);

        double res = cumulValues.get(e).getValue();
        res -= cumulValues.get(b - 1).getValue();
        return res;
    }

    /**
     * Gets the Squared Mean for elements within index b..e, multiplied by weight. Note that
     * n*mean^2 = sum^2/n when mean := sum/n
     *
     * @param b index of begin element
     * @param e index of end element
     * @return the sum of squared mean
     */
    private double getSSM(int b, int e) {
        double res = getSumOfWeightedValues(b, e);
        return res * res / getSumOfWeights(b, e);
    }

    /**
     * Finds CB[i+completedRows] given that the result is at least
     * bp+(completedRows-1) and less than ep+(completedRows-1)
     * Complexity: O(ep-bp) <= O(m) @
     *
     * @param i  startIndex
     * @param bp endindex
     * @param ep
     * @return the index
     */
    private int findMaxBreakIndex(int i, int bp, int ep) {
        assert (bp < ep);
        assert (bp <= i);
        assert (ep <= i + 1);
        assert (i < bufferSize);
        assert (ep <= bufferSize);

        double minSSM = previousSSM[bp] + getSSM(bp + completedRows, i + completedRows);
        int foundP = bp;
        while (++bp < ep) {
            double currSSM = previousSSM[bp] + getSSM(bp + completedRows, i + completedRows);
            if (currSSM > minSSM) {
                minSSM = currSSM;

                foundP = bp;
            }
        }
        currentSSM[i] = minSSM;
        return foundP;
    }

    /**
     * Find CB[i+completedRows] for all i>=bi and i<ei given that the
     * results are at least bp+(completedRows-1) and less than
     * ep+(completedRows-1)
     * Complexity: O(log(ei-bi)*Max((ei-bi),(ep-bp)))
     * <= O(m*log(m))
     *
     * @param bi
     * @param ei
     * @param bp
     * @param ep
     * @return
     */
    private void calcRange(int bi, int ei, int bp, int ep) {
        assert (bi <= ei);

        assert (ep <= ei);
        assert (bp <= bi);

        if (bi == ei) {
            return;
        }
        assert (bp < ep);

        int mi = (bi + ei) / 2;
        int mp = findMaxBreakIndex(mi, bp, Math.min(ep, mi + 1));

        assert (bp <= mp);
        assert (mp < ep);
        assert (mp <= mi);

        // solve first half of the sub-problems with lower 'half' of possible outcomes
        calcRange(bi, mi, bp, Math.min(mi, mp + 1));
        classBreaks[classBreaksIndex + mi] = mp; // store result for the middle element.
        // solve second half of the sub-problems with upper 'half' of possible outcomes
        calcRange(mi + 1, ei, mp, ep);
    }

    /**
     * Swaps the content of the two lists with each other.
     */
    private void swapArrays() {
        double[] temp = previousSSM;
        previousSSM = currentSSM;
        currentSSM = temp;
    }

    /**
     * Starting point of calculation of breaks.
     * <p>
     * complexity: O(m*log(m)*k)
     */
    private void calcAll() {
        if (numBreaks >= 2) {
            classBreaksIndex = 0;
            for (completedRows = 1; completedRows < numBreaks - 1; ++completedRows) {
                calcRange(0, bufferSize, 0, bufferSize); // complexity: O(m*log(m))
                swapArrays();
                classBreaksIndex += bufferSize;
            }
        }
    }

    @RequiredArgsConstructor
    @Data
    private static class ClassifierAndBreaks {
        private final int[] breaksIndices;
        private final JenksFisher jf;
    }

    /**
     * Does the internal processing to actually create the breaks.
     */
    public static double[] classifyJenksFisherFromValueCountPairs(int k, List<ValueCountPair> vcpc) {
        return Arrays.stream(calculateBreaks(k, vcpc).breaksIndices)
                .mapToDouble(i -> vcpc.get(i).getValue())
                .toArray();
    }

    private static ClassifierAndBreaks calculateBreaks(int k, List<ValueCountPair> vcpc) {
        int[] breaksArray = new int[k];
        int m = vcpc.size();
        assert (k <= m); // PRECONDITION
        if (k == 0)
            return new ClassifierAndBreaks(breaksArray, null);
        final JenksFisher jf = new JenksFisher(k, vcpc);
        if (k > 1) {
            // runs the actual calculation
            jf.calcAll();
            int lastClassBreakIndex = jf.findMaxBreakIndex(jf.bufferSize - 1, 0, jf.bufferSize);
            while (--k != 0) {
                // assign the break values to the result
                breaksArray[k] = lastClassBreakIndex + k;
                assert (lastClassBreakIndex < jf.bufferSize);
                if (k > 1) {
                    jf.classBreaksIndex -= jf.bufferSize;
                    lastClassBreakIndex = jf.classBreaks[jf.classBreaksIndex + lastClassBreakIndex];
                }
            }
            assert (jf.classBreaks[jf.classBreaksIndex] == jf.classBreaks[0]);
        }
        assert (k == 0);
        breaksArray[0] = 0; // break for the first class is the minimum of the dataset.
        return new ClassifierAndBreaks(breaksArray, jf);
    }

    /**
     * Main entry point for creation of Jenks-Fisher natural breaks.
     *
     * @param values array of the values, do not need to be sorted.
     * @param k      number of breaks to create
     * @return Array with breaks
     */
    static List<Double> createJenksFisherBreaksArray(List<Double> values, int k) {
        List<ValueCountPair> sortedUniqueValueCounts = getValueCountPairs(values);
        return sortedUniqueValueCounts.size() > k ?
                Arrays.stream(classifyJenksFisherFromValueCountPairs(k, sortedUniqueValueCounts))
                        .boxed()
                        .collect(Collectors.toList()) :
                sortedUniqueValueCounts.stream().mapToDouble(ValueCountPair::getValue)
                        .boxed()
                        .collect(Collectors.toList());
    }

    /**
     * Calculates the occurrence count of given values and returns them.
     *
     * @param values
     * @return Occurences of values.
     */
    public static List<ValueCountPair> getValueCountPairs(List<Double> values) {
        return values
                .stream()
                .collect(Collectors.groupingBy(v->v, Collectors.counting()))
                .entrySet()
                .stream()
                .map(e->new ValueCountPair(e.getKey(), e.getValue()))
                .sorted()
                .collect(Collectors.toList());
    }

    /**
     * Returns an ordered list of clusters based on break values, each with a pre-calculated variance value, size(items per cluster) and center
     *
     * @param k    number of breaks
     * @param vcpc asc ordered input of values and their occurence counts.
     * @return the ordered clusters
     */
    public static List<Cluster> classifyJenksFisherClustersFromValueCountPairs(final int k, List<ValueCountPair> vcpc) {
        ClassifierAndBreaks classifierAndBreaks = calculateBreaks(k, vcpc);
        final JenksFisher jf = classifierAndBreaks.jf;
        List<Cluster> result = new ArrayList<>(k);
        int[] breakIndices = classifierAndBreaks.breaksIndices;
        ValueCountPair lower = null;
        for (int i = 0; i < k; ++i) {
            ValueCountPair upper = jf.cumulValues.get(i == k - 1 ? jf.cumulValues.size() - 1 : breakIndices[i + 1] - 1);
            double deltaValue = lower == null ? upper.getValue() : upper.getValue() - lower.getValue();
            long deltaCount = lower == null ? upper.getCount() : upper.getCount() - lower.getCount();
            double center = deltaValue / deltaCount;

            // Calculate variance
            final int vcpcStartIndex = breakIndices[i];
            final int vcpcEndIndex = (i == k - 1) ? vcpc.size() : breakIndices[i + 1];
            List<ValueCountPair> subList = vcpc.subList(vcpcStartIndex, vcpcEndIndex);
            double variance = new Variance().evaluate(subList.stream()
                            .mapToDouble(v -> Math.abs(v.getValue() - center)).toArray(),
                    subList.stream().mapToDouble(ValueCountPair::getCount).toArray());
            result.add(new Cluster(center, variance, deltaCount));
            lower = upper;
        }
        return result;
    }
}
