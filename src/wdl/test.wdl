version 1.0

    #Task definition

    task say_hello {
        input {
            String name    
        }
        command {
        set -euxo pipefail
        echo "hello ~{name}"
        echo "hello ~{name} > greetings.txt"
        }
        
        output{
            File greeting = "greetings.txt"
            
        }
        runtime {
            docker: "debian:stretch-slim"    
        }
    }

    workflow hello {
        input {
            String name    
        }    
        call say_hello {
            input:
                name = name
        }

        output {
            File greeting = say_hello.greeting    
        }
    }

