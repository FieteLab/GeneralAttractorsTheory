module ProjectSupervisor
    export Supervisor, 
        move_to_datadir, 
        store_data,
        set_datadir,
        get_entries,
        fetch,
        save_plot,
        generate_or_load

    using ObjectivePaths, DataFrames, Term, MyterialColors, FileIO, UUIDs, Dates, LibGit2, YAML
    using Term.Tables
    import Plots: savefig

    # ---------------------------------------------------------------------------- #
    #                                  SUPERVISOR                                  #
    # ---------------------------------------------------------------------------- #

    # --------------------------------- creation --------------------------------- #
    """
        Supervisor

    Keeps an eye on your project. Storing metadata on what gets saved where.
    """
    mutable struct Supervisor
        projectdir::Folder
        datadir::Folder
        gitrepo::Union{Nothing, GitRepo}
        metadatafile::ObjectivePaths.File
        metadata::AbstractDict
    end

    function Supervisor(project_name::String)
        fld = Folder(pwd()) / project_name
        while !isnothing(fld) && !exists(fld)
            fld = (fld - 2) / project_name
        end
        isnothing(fld) && error("Could not find project folder $project_name")

        # get/create a data dir
        datadir = fld / "data"
        exists(datadir) || mkdir(datadir)

        meta, meta_path = load_or_create_metadata(datadir)
        Supervisor(fld, datadir, get_gitrepo(fld), meta_path, meta)
    end


    # repr for Supervisor showing folder
    function Base.show(io::IO, supervisor::Supervisor)
        meta = update_metadata!(supervisor.metadata)
        tprintln(io, "{bold $orange}Supervisor with datadir:{/bold $orange}")
        println(io, supervisor.datadir)

        tprintln(io, "{$orange}Git repository:{/$orange}")
        if isnothing(supervisor.gitrepo)
            println(io, "No git repository found")
        else
            println(io, supervisor.gitrepo)
        end

        df = DataFrame(meta)
        tprintln(io, "{$orange}Metadata $(size(df)):{/$orange}"; highlight=false)
        "tag" ∉ keys(meta) && return
        _meta = filter(k -> k.first ∈ ("date","time","folder", "tag"), meta)

        # keep at most 5 entries
        for k in keys(_meta)
            length(_meta[k]) < 5 && continue
            _meta[k] = rand(_meta[k], 5)
        end

        println()
        tprintln(io, Table(_meta;
            columns_style=["bold white", "default", "dim", "dim"],
            footer = length,
            footer_style = "blue",
            header_style = "blue bold",
        ))

        println.(keys(supervisor.metadata))
    end
    

    """
    Find a folder that is a Git repo.
    Either the project or one of its parents.
    """
    function get_gitrepo(fld::Folder)::Union{Nothing, GitRepo}
        repo = nothing
        repodir = fld
        while isnothing(repo)
            repo = try
                LibGit2.GitRepoExt(repodir.path)
            catch e
                nothing
            end
            repodir = repodir - 1
            isnothing(repodir) && break
            exists(repodir) || break
        end
        return repo
    end

    """
        Get info for current git repository (e.g. commit)
    """
    function git_info(sup::Supervisor)
        repo = sup.gitrepo
        isnothing(repo) && return "No git repository found"

        suffix = LibGit2.isdirty(repo) ? "-dirty" : ""
        c = try
            gdr = LibGit2.GitDescribeResult(repo)
            fopt = LibGit2.DescribeFormatOptions(dirty_suffix=pointer(suffix))
            LibGit2.format(gdr, options=fopt)
        catch GitError
            string(LibGit2.head_oid(repo)) * suffix
        end
        return c
    end



    # ------------------------------- change folder ------------------------------ #
    
    """
        Move datadir to a subfolder
    """
    function move_to_datadir(supervisor::Supervisor, fld::String)
        fld = supervisor.datadir / fld
        exists(fld) || mkdir(fld)
        supervisor.datadir = fld
        meta, meta_path = load_or_create_metadata(fld)
        supervisor.metadata = meta
        supervisor.metadatafile = meta_path
    end

    """
        Move datadir to another folder.
    """
    function set_datadir(supervisor::Supervisor, fld::String)
        fld = Folder(fld)
        supervisor.datadir = fld
        meta, meta_path = load_or_create_metadata(fld)
        supervisor.metadata = meta
        supervisor.metadatafile = meta_path
    end


    # --------------------------------- metadata --------------------------------- #

    function validate(meta::AbstractDict)
        @assert eltype(keys(meta)) == String "Metadata keys must be Strings: $(eltype(keys(meta)))"
        # make sure all entries have the same length
        if "uid" ∈ keys(meta)
            n_entries = unique(length.(values(meta)))
            _counts = map(p -> (p.first, length(p.second)), collect(pairs(meta)))
            @assert length(n_entries) <= 1 "Metadata entries have different lengths $n_entries- $_counts"
        end
    end

    function load_or_create_metadata(fld::Folder)::Tuple{AbstractDict, ObjectivePaths.File}
        path = fld / "metadata.yaml"
        meta = if exists(path)
            YAML.load_file(path.path) |> update_metadata!
        else
            Dict()
        end
        meta = Dict{String, Any}(meta)

        # make sure vecs can accept any type
        for (k, v) in pairs(meta)
            length(v) == 0 && continue
            meta[string(k)] = Any[v...]
        end

        validate(meta)
        return meta, path
    end

    """
        update_metadata!(metadata::AbstractDict)::AbstractDict

    Remove obsolete entries pointing to files that no longer exist.
    """
    function update_metadata!(metadata::AbstractDict)::AbstractDict
        "name" ∉ keys(metadata) && return metadata
        to_remove = []
        for (i, f) in enumerate(metadata["name"])
            path = Folder(metadata["folder"][i]) / f
            exists(path) && continue
            
            push!(to_remove, i)
        end

        for k in keys(metadata)
            metadata[k] = deleteat!(metadata[k], to_remove)
        end

        length(metadata["name"]) == 0 && (metadata = Dict{String, Any}())
        metadata = Dict{String, Any}(metadata)
        validate(metadata)
        return metadata
    end

    # """

    # Check if any entry in supervisor.metadata matches the given metadata.
    # If so, return the uid of the entry, otherwise return nothing
    # """
    # function check_entry_exists(sup::Supervisor, metadata::AbstractDict, folder::Folder)::Union{Nothing, String}
    #     uid = ""
    #     for (k,v) in pairs(metadata)
    #         k = string(k)

    #         (v isa AbstractArray &&
    #             k ∉ keys(sup.metadata) &&
    #             v ∉ sup.metadata[k] &&
    #             length(sup.metadata[k]) == 0 ) && return nothing

    #         matches = sup.metadata[k] .== v
    #         any(Bool.(matches)) || return nothing

    #         fld = sup.metadata["folder"][sup.metadata[k] .== v][1]
    #         fld == folder.path || return nothing
    #         uid = sup.metadata["uid"][sup.metadata[k] .== v][1]
    #     end
    #     return uid
    # end

    function add_metadata_entry(
        sup::Supervisor,
        uid::String,
        date::SubString,
        time::SubString,
        git_info::String,
        savepath::ObjectivePaths.File,
        metadata::AbstractDict
    )
        entry = Dict(
            "uid" => uid,
            "date" => date,
            "time" => time,
            "git_info" => git_info,
            "folder" => dirname(savepath),
            "name" => name(savepath),
            pairs(metadata)...
        )

        n_prev_entries =  haskey(sup.metadata, "uid") ? length(sup.metadata["uid"]) : 0
        for (k, v) in pairs(entry)
            k = string(k)

            # if metadata doesn't have this key, add it to the metadata
            if k ∉ keys(sup.metadata)
                sup.metadata[k] = repeat(Any[""], n_prev_entries)
            end
            
            # store value
            push!(sup.metadata[k], v)
        end

        # make sure we fill metadata keys not in this entry
        ekeys = string.(keys(entry))
        for k in string.(keys(sup.metadata))
            k ∈ ekeys && continue
            push!(sup.metadata[k], "")
        end

        validate(sup.metadata)
    end

    function get_entries(sup::Supervisor; filter_values...)
        df = DataFrame(sup.metadata)
        for (k, v) in pairs(filter_values)
            k = string(k)
            k ∉ keys(sup.metadata) && continue
            df = df[df[:, k] .== v, :]
        end
        return df
    end

    # --------------------------------- save data -------------------------------- #

    """
        get_savename(sup, name, fmt, paths_elements...; )

    Get the path where some data will be saved.
    """
    function get_savename(sup, name, fmt, paths_elements...; )
        # get path to destination folder
        dest = if length(paths_elements) == 0
            sup.datadir
        else
            folder = joinpath(paths_elements...)
            _dest = sup.datadir / folder
            exists(_dest) || mkdir(_dest)
            _dest
        end

        # get save name
        name = string(name) * "." * fmt
        return dest/name
    end
        
    """
        generate_or_load(
            fn, 
            sup::Supervisor,
            paths_elements...;
            name=nothing,
            fmt=nothing,
            kwargs...
            )

        If the file exists, load it. Otherwise, generate it using `fn` and save it.
    """
    function generate_or_load(
        fn, 
        sup::Supervisor,
        paths_elements...;
        name=nothing,
        fmt=nothing,
        load_existing = true,
        kwargs...
        )

        savename = get_savename(sup, name, fmt, paths_elements...)
        if exists(savename)
            return if load_existing
                load(savename.path)
            else
                nothing
            end
        else
            data = fn()
            store_data(sup; savename=savename, data=data, kwargs...)
        end
    end

    """
        store_data(
            sup::Supervisor,
            data;
            paths_elements...;
            name=nothing,
            fmt=nothing,
            kwargs...
        )

    Store some data in the supervisor's data folder.
    """
    function store_data(
        sup::Supervisor,
        paths_elements...;
        data = nothing,
        name=nothing,
        fmt=nothing,
        kwargs...
    )
        savename = get_savename(sup, name, fmt, paths_elements...)
        return store_data(sup; savename=savename, data=data, kwargs...)
    end

    """
        store_data(
            sup::Supervisor,
            savename,
            data;
            metadata::AbstractDict = Dict("tag"=>"notag"),
        )

    Store some data in the supervisor's data folder.
    """
    function store_data(
        sup::Supervisor;
        savename::ObjectivePaths.File,
        data,
        metadata::AbstractDict = Dict("tag"=>"notag"),
    )
        @assert "tag" ∈ string.(keys(metadata)) "Data entry metadata must have a tag."

        # collect meta info
        uid = string(uuid4())
        when = Dates.now()
        date, time = split(string(when), "T")


        try
            save(savename.path, data)
        catch e
            @error "Error saving data at $(savename.path)" e
        end

        # add metadata entry
        add_metadata_entry(sup, uid, date, time, git_info(sup), savename, metadata)

        # save metadata
        YAML.write_file(sup.metadatafile.path, sup.metadata)
    end


    function save_plot(sup::Supervisor, plt, name::AbstractString; dpi=300)
        meta = git_info(sup)
        path = (sup.projectdir / "plots" / (name*"_"*meta)).path
        
        savefig(plt, path * ".png")
        savefig(plt, path * ".svg")
    end
    


    # --------------------------------- load data -------------------------------- #

    """
        fetch(sup::Supervisor; filter_criteria...)

    Fetch a subset of the data stored in the metadata based 
    on some criteria to filter which entries to select.
    """
    function fetch(sup::Supervisor; filter_criteria...)::Tuple{DataFrame, Vector}
        df = get_entries(sup; filter_criteria...)
        
        data = []
        for entry in eachrow(df)
            path = joinpath(entry.folder, entry.name)
            push!(data, load(path))
        end

        return df, data
    end
end