package org.rcsb.biozernike.ligzernike;
import picocli.CommandLine;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class DoubleListConverter implements CommandLine.ITypeConverter<List<Double>> {
    @Override
    public List<Double> convert(String value) {
        return Arrays.stream(value.split(","))
                .map(Double::parseDouble)
                .collect(Collectors.toList());
    }
}