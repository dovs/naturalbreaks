package de.pschoepf.naturalbreaks;

import lombok.Data;
import lombok.RequiredArgsConstructor;
import lombok.ToString;

@RequiredArgsConstructor
@Data
@ToString
public class Cluster {
    private final double center;
    private final double variance;
    private final long size;
}
