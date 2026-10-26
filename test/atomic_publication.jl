# Keep the source readable but not renameable after write(path, payload) closes it.
struct LockedCheckpointPayload
    directory::String
    handle::Base.RefValue{Ptr{Cvoid}}
end
function Base.write(io::IO, payload::LockedCheckpointPayload)
    n = write(io, b"checkpoint")
    flush(io)
    source = only(filter(p -> endswith(p, ".partial"), readdir(payload.directory; join = true)))
    handle = ccall((:CreateFileW, "kernel32"), stdcall, Ptr{Cvoid},
        (Cwstring, UInt32, UInt32, Ptr{Cvoid}, UInt32, UInt32, Ptr{Cvoid}),
        source, 0x80000000, 0x00000003, C_NULL, 0x00000003, 0x00000000, C_NULL)
    handle == Ptr{Cvoid}(typemax(UInt)) && error("could not lock the temporary checkpoint")
    payload.handle[] = handle
    n
end

@testset "Checkpoint rename never falls back to copying" begin
    mktempdir() do dir
        destination = joinpath(dir, "checkpoint.jsonl")
        if Sys.iswindows()
            payload = LockedCheckpointPayload(dir, Ref(Ptr{Cvoid}(C_NULL)))
            try
                @test_throws Base.IOError BOLDassignment._save_atomic(destination, payload)
                @test !isfile(destination)
            finally
                if payload.handle[] != C_NULL
                    ccall((:CloseHandle, "kernel32"), stdcall, Cint,
                          (Ptr{Cvoid},), payload.handle[])
                end
            end
        end
        BOLDassignment._save_atomic(destination, b"first")
        @test read(destination, String) == "first"
        BOLDassignment._save_atomic(destination, b"second")
        @test read(destination, String) == "second"
    end
end
