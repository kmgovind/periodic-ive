module Notify

using HTTP

export send_ntfy

function send_ntfy(message, title="ASV Simulation Update", priority="default")
    # Replace 'my_secret_asv_sim_topic' with a unique string of your choice
    topic = "my_secret_asv_sim_topic" 
    try
        HTTP.put("https://ntfy.sh/$topic", 
                 headers=[
                     "Title" => title, 
                     "Priority" => priority,
                     "Tags" => priority == "high" ? "warning,skull" : "checkered_flag,ocean"
                 ],
                 body=message)
    catch e
        println("Notification failed: $e")
    end
end

end # End of Notify module