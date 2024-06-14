rule tab2xlsx:
    input: 
        "{file}"
    output: 
        "{file}.xlsx"
    shell: 
        "tab2xlsx < {input} > {output}"