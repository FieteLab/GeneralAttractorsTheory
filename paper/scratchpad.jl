using YAML

source = "F:\\PostDoc\\GeneralAttractors\\data\\mfld_top2\\metadata.yaml"
old = "/Users/federicoclaudi/Desktop/GeneralAttractors/data/mfld_top2/"
new = "F:\\PostDoc\\GeneralAttractors\\data\\mfld_top2\\"


function clean_yaml()
    # load adata
    data = YAML.load_file(source)

    # for each entry, change the folder 
    for (key, value) in data
        # get the last name of the foplder path
        data[key]["folder"] = replace(value["folder"], old => new)
    end

    # save the new data
    YAML.write_file("F:\\PostDoc\\GeneralAttractors\\data\\mfld_top2\\metadata.yaml", data)

end

clean_yaml()