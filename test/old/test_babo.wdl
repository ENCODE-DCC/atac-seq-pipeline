task T {
    String? s
    command {
        echo ${s}
    }
    output {
        String ret = stdout()
    }
}

workflow X {
    Boolean? flag

    call T { input: s = if select_first([flag,false]) then "OKAY" else "FAIL"}
}
