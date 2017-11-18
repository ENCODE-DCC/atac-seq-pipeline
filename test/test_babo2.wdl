workflow X {
    Boolean? flag
    output {
        String s = if select_first([flag,false]) then "OKAY" else "FAIL"
    }
}
