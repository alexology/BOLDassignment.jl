# Times of the last requests. Shared by every engine built with the default
# `throttle`, because the BOLD limits apply to the client, not to one call.
mutable struct _Throttle
    lock::ReentrantLock
    last_submit::Float64
    last_request::Float64
end
_Throttle() = _Throttle(ReentrantLock(), -Inf, -Inf)

const _SHARED_THROTTLE = _Throttle()

"""
    IDEngine(; kwargs...)

Connection settings for the BOLD Barcode ID API.

The API documents at most 3 submissions and 40 status or result requests per
minute; going over may block the client for a while. The defaults stay below
both limits. The spacing is shared by all engines using the default
`throttle`, so consecutive or concurrent calls respect the limits together.

# Keywords
- `base_url = "https://id.boldsystems.org"`
- `submit_interval = 21.0`: minimum seconds between two submissions.
- `request_interval = 1.6`: minimum seconds between two status/result requests.
- `status_wait = 15.0`: seconds between two status checks of the same submission.
- `max_active = 5`: submissions waiting for results at the same time.
- `timeout = 3600.0`: seconds after which an unfinished submission is resubmitted.
- `max_resubmits = 3`: resubmissions of one chunk before giving up.
- `max_retries = 8`: consecutive failures tolerated for one request.
- `retry_wait = 30.0`: wait after a failure, multiplied by the attempt number.
- `request`: function `(method, url, payload) -> (status, body)` performing the
  HTTP request; `payload` is `nothing` or the FASTA bytes. Replaceable in tests.
- `sleep`, `time`: clock functions, replaceable in tests.
- `throttle`: request spacing state; pass `BOLDassignment._Throttle()` for an
  independent one (tests).
- `verbose = true`: print progress messages.
"""
Base.@kwdef struct IDEngine
    base_url::String = "https://id.boldsystems.org"
    submit_interval::Float64 = 21.0
    request_interval::Float64 = 1.6
    status_wait::Float64 = 15.0
    max_active::Int = 5
    timeout::Float64 = 3600.0
    max_resubmits::Int = 3
    max_retries::Int = 8
    retry_wait::Float64 = 30.0
    request::Function = _http_request
    sleep::Function = Base.sleep
    time::Function = Base.time
    throttle::_Throttle = _SHARED_THROTTLE
    verbose::Bool = true
    function IDEngine(base_url, submit_interval, request_interval, status_wait, max_active,
                      timeout, max_resubmits, max_retries, retry_wait, request, sleep, time,
                      throttle, verbose)
        for (name, v) in (("submit_interval", submit_interval), ("request_interval", request_interval),
                          ("status_wait", status_wait), ("timeout", timeout),
                          ("retry_wait", retry_wait))
            (isfinite(v) && v >= 0) || throw(ArgumentError("$name must be finite and non-negative"))
        end
        for (name, v) in (("max_active", max_active), ("max_retries", max_retries))
            v >= 1 || throw(ArgumentError("$name must be at least 1"))
        end
        max_resubmits >= 0 || throw(ArgumentError("max_resubmits must be non-negative"))
        new(rstrip(base_url, '/'), submit_interval, request_interval, status_wait, max_active,
            timeout, max_resubmits, max_retries, retry_wait, request, sleep, time, throttle,
            verbose)
    end
end

function _http_request(method::String, url::String, payload)
    headers = ["User-Agent" => "BOLDassignment.jl"]
    response = if payload === nothing
        HTTP.request(method, url, headers; status_exception = false, retry = false,
                     readtimeout = 600)
    else
        form = HTTP.Form(Dict("fasta_file" => HTTP.Multipart("query.fasta", IOBuffer(payload),
                                                             "text/plain")))
        HTTP.request(method, url, headers, form; status_exception = false, retry = false,
                     readtimeout = 600)
    end
    return (Int(response.status), response.body)
end

function _wait_until(engine::IDEngine, t::Float64)
    dt = t - engine.time()
    dt > 0 && engine.sleep(dt)
    return nothing
end

_log(engine::IDEngine, msg) = engine.verbose && @info msg

# One request with spacing and retries. `accept` lists statuses returned to
# the caller instead of being retried (e.g. 404 for unfinished results).
function _call(engine::IDEngine, method::String, url::String, payload = nothing;
               accept = (200,))
    submit = payload !== nothing
    throttle = engine.throttle
    for attempt in 1:engine.max_retries
        status, body = lock(throttle.lock) do
            if submit
                _wait_until(engine, throttle.last_submit + engine.submit_interval)
            end
            _wait_until(engine, throttle.last_request + engine.request_interval)
            try
                engine.request(method, url, payload)
            catch err
                err isa InterruptException && rethrow()
                (-1, Vector{UInt8}(codeunits(sprint(showerror, err))))
            finally
                throttle.last_request = engine.time()
                submit && (throttle.last_submit = throttle.last_request)
            end
        end
        status in accept && return status, body
        if 400 <= status < 500 && status != 429 && status != 408
            throw(ErrorException("BOLD request failed with HTTP $status: $url\n" *
                                 _snippet(body)))
        end
        _log(engine, "BOLD request failed (HTTP $status, attempt $attempt/$(engine.max_retries)); retrying")
        attempt < engine.max_retries && engine.sleep(engine.retry_wait * attempt)
    end
    throw(ErrorException("BOLD request failed $(engine.max_retries) times: $url"))
end

_snippet(body) = first(String(copy(body)), 500)

_submission_url(engine::IDEngine, db::String, mode) =
    string(engine.base_url, "/submission?db=", HTTP.escapeuri(db), "&mi=", mode.mi,
           "&mo=", mode.mo, "&maxh=", mode.maxh)

function _submit(engine::IDEngine, db::String, mode, fasta::Vector{UInt8})
    url = _submission_url(engine, db, mode)
    for attempt in 1:engine.max_retries
        _, body = _call(engine, "POST", url, fasta)
        obj = try
            JSON3.read(body)
        catch
            nothing
        end
        if obj isa JSON3.Object && get(obj, :sub_id, nothing) isa AbstractString
            return String(obj[:sub_id])
        end
        # The engine answers with a `detail` field when the usage limit is reached.
        _log(engine, "BOLD submission refused (attempt $attempt): $(_snippet(body))")
        attempt < engine.max_retries && engine.sleep(engine.retry_wait * attempt)
    end
    throw(ErrorException("BOLD refused the submission $(engine.max_retries) times"))
end

# Whether the submission is finished. While running the engine answers with
# counts (`{"queued":0,"processing":2,"completed":0}`), when finished with
# `{"status":"complete"}` (observed 2026-10-23). Unknown or unreadable
# answers count as not finished; the results are validated anyway.
function _finished(engine::IDEngine, sub_id::String)
    url = string(engine.base_url, "/submission/status/", HTTP.escapeuri(sub_id))
    status, body = _call(engine, "GET", url; accept = (200, 404))
    status == 404 && return false
    obj = try
        JSON3.read(body)
    catch
        nothing
    end
    obj isa JSON3.Object || return false
    st = get(obj, :status, nothing)
    st isa AbstractString && return lowercase(st) in ("complete", "completed")
    count(k) = (v = get(obj, k, nothing); v isa Integer ? Int(v) : -1)
    queued, processing, completed = count(:queued), count(:processing), count(:completed)
    return queued == 0 && processing == 0 && completed > 0
end

# Raw JSON Lines body, or `nothing` if not ready yet.
function _results(engine::IDEngine, sub_id::String)
    url = string(engine.base_url, "/submission/results/", HTTP.escapeuri(sub_id))
    status, body = _call(engine, "GET", url; accept = (200, 404))
    return status == 404 ? nothing : body
end
