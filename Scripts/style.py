import streamlit as st


def style():
    """
    Apply custom styling to the Streamlit app, including page configuration and CSS modifications.
    """

    # Hide the hamburger menu and footer
    hide_st_style = """
        <style>
            #Mainmenu {visibility: hidden;}
            header{visibility: visible;}
        <style>
    """
    st.markdown(hide_st_style, unsafe_allow_html=True)

    # Custom CSS code
    custom_style = """
        <style>
            ::-webkit-scrollbar {
                background: rgb(0, 170, 240);
                height: 6px;
                width: 5px;
            }

            body * {
                font-family: 'Arial', sans-serif !important;
            }

            /* Ensure the sidebar is an overlay and does not affect the main content */
            section[data-testid="stSidebar"] {
                position: fixed !important; /* Fix the sidebar in place */
                width: 325px !important; /* Fixed width */
                height: 100% !important; /* Full height */
                left: 0 !important; /* Align to the left */
                transition: transform 0.3s ease-in-out; /* Smooth transition when toggling */
            }

            /* Handle sidebar hidden state */
            section[data-testid="stSidebar"][aria-expanded="false"] {
                transform: translateX(-100%) !important; /* Hide the sidebar off-screen */
            }

            /* Main content adjustments */
            section[data-testid="stMain"] {
                margin-left: 325px !important; /* Reserve space for the sidebar */
                transition: margin-left 0.3s ease-in-out; /* Smooth adjustment when toggling */
            }

            /* When sidebar is hidden */
            section[data-testid="stSidebar"][aria-expanded="false"] ~ section[data-testid="stMain"] {
                margin-left: 0 !important; /* Adjust content to fill the space */
            }

            /* Ensure no extra padding/margin around the content */
            .block-container, div[data-testid="stAppViewContainer"] {
                padding-top: 0rem !important;
                padding-bottom: 0rem !important;
                margin-top: 0rem !important;
                margin-left: 0rem !important; /* Align closer to the sidebar */
            }

            /* Reduce spacing between the sidebar and content */
            section[data-testid="stSidebar"] {
                margin-right: -4rem !important; /* Ensure no space on the right of the sidebar */
            }


            /* Target the FileUploader drop zone */
            section[data-testid="stFileUploaderDropzone"] {
                border-color: rgb(0, 0, 0) !important;
                padding: 0.25rem !important;
                background-color: #00aaff57 !important; /* Light blue background */
                border-radius: 0.50rem !important;
                color: rgb(0, 0, 0) !important;
                border: 1.5px solid rgb(0, 0, 0) !important;
                display: flex !important; /* Use flexbox */
                flex-direction: column !important;
                justify-content: center !important; /* Center vertically */
                align-items: center !important; /* Center horizontally */
                height: 100px;
            }

            /* Target the FileUploader drop zone button */
            section[data-testid="stFileUploaderDropzone"] button {
                width: 250px !important; /* Set desired width */
                margin-top: -15px !important; /* Reset margins for centering */
            }

            /* Target the FileUploader drop zone button */
            div[data-testid="stFileUploaderFile"] {
                margin-top: -9px !important; /* Reset margins for centering */
            }


            /* General button style */
            button[data-testid="stBaseButton-secondary"],
            button[data-testid="stBaseButton-primary"],
            button[data-testid="stBaseButton-secondaryFormSubmit"],
            button[data-testid="stBaseButton-upload"],
            button[data-testid="stBaseButton-download"] {
                display: inline-flex !important;
                -webkit-box-align: center !important;
                align-items: center !important;
                -webkit-box-pack: center !important;
                justify-content: center !important;
                font-weight: 500 !important;
                padding: 0.5rem 1rem !important; /* Slightly larger padding */
                border-radius: 0.5rem !important;
                min-height: 2.5rem !important;
                margin: 0px !important;
                line-height: 1.6 !important;
                color: white !important; /* White text for contrast */
                width: auto !important;
                user-select: none !important;
                background-color: rgb(0, 153, 255) !important; /* Darker blue */
                border: 2px solid rgb(255 0 0) !important; /* Light blue border */
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2); /* Subtle shadow */
                transition: all 0.3s ease-in-out; /* Smooth hover effect */
            }


            /* Hover effect */
            button[data-testid="stBaseButton-secondary"]:hover,
            button[data-testid="stBaseButton-primary"]:hover,
            button[data-testid="stBaseButton-secondaryFormSubmit"]:hover,
            button[data-testid="stBaseButton-upload"]:hover,
            button[data-testid="stBaseButton-download"]:hover {
                background-color: rgb(0, 102, 204) !important; /* Slightly lighter blue */
                border-color: rgb(163 0 0) !important; /* Softer light blue */
                box-shadow: 0 5px 8px rgba(0, 0, 0, 0.4); /* Slightly stronger shadow */
                transform: translateY(-2px); /* Lift effect */
            }

                /* Style for disabled buttons */
            button[data-testid="stBaseButton-secondary"][disabled] {
                background-color: rgb(200, 200, 200) !important; /* Static gray background */
                color: rgb(120, 120, 120) !important; /* Lighter gray text */
                border: 2px solid rgb(180, 180, 180) !important; /* Slightly darker border */
                box-shadow: none !important; /* Remove shadow */
                cursor: not-allowed !important; /* Indicate the button is not clickable */
                transform: none !important; /* No movement on hover */
                transition: none !important; /* No transition effects */
            }

            [data-testid="stAlert"] {
                border: 1.5px solid rgb(255, 0, 0);
                border-radius: 15px;
                font-family: Arial, sans-serif;
                font-size: 16px;
                color: black;
                box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.35);
                transition: transform 0.5s ease, box-shadow 0.5s ease, background-color 0.5s ease;
                background-color: rgba(255, 255, 255, 0.8);
            }

            /* Hover Effect */
            [data-testid="stAlert"]:hover {
                transform: scale(1.02);
                box-shadow: 0px 6px 12px rgba(255, 0, 0, 0.6);
                background-color: rgba(255, 245, 245, 0.9);
            }

            /* Subtle Pulse Animation */
            [data-testid="stAlert"] {
                animation: pulse 2s infinite;
            }

            @keyframes pulse {
                0% {
                    box-shadow: 0px 2px 4px rgba(0, 0, 0, 0.2);
                }
                50% {
                    box-shadow: 0px 3px 5px rgba(255, 0, 0, 0.3);
                }
                100% {
                    box-shadow: 0px 2px 4px rgba(0, 0, 0, 0.2);
                }
            }

            /* Target the specific question mark (help icon) */
            div[data-testid="stTooltipHoverTarget"] {
                transform: none !important; /* Prevent size changes for non-hover state */
                transition: all 0.3s ease-in-out; /* Smooth transition effect */
                display: inline-flex !important;
                align-items: center !important;
                justify-content: center !important;
            }

            /* Style the question mark icon (red color) */
            div[data-testid="stTooltipHoverTarget"] svg {
                stroke: red !important; /* Default red color */
                transition: stroke 0.3s ease-in-out; /* Smooth color transition */
                z-index: 2 !important;
            }

            /* Hover effect for the question mark */
            div[data-testid="stTooltipHoverTarget"]:hover svg{
                transform: scale(1.2) !important; /* Slightly increase size on hover */
                box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2) !important; /* Add shadow */
            }

            div[data-testid="stToast"] {
                margin-top: 37px;
                border: 2px solid rgb(255, 0, 0);
                border-image: linear-gradient(to right, rgb(255, 0, 0), rgb(0, 0, 255)) 1;
                animation: pulse 2s infinite;
                transform: scale(1.075)    
                font-size: 16px;  
            }

            div[data-testid="stTooltipHoverTarget"]:hover svg {
                stroke: darkred !important; /* Change stroke color on hover */
            }


            /* Style the tab container for sticky behavior */
            div[data-baseweb="tab-list"] {
                display: flex;
                position: sticky;
                table-layout: auto;
                justify-content: space-between;; /* Align tabs to the left */
                gap: 0.1px; /* Smaller gap between tabs */
                background-color: #ffffff; /* Light blue background matching your theme */
                padding: 0.1rem 0.1rem !important;
                border-radius: 0.4rem; /* Slightly rounded edges */
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0); /* Subtle shadow for depth */
                margin-top: -5% !important; /* Pull the tabs higher */
                margin-left: auto; /* Center tabs in container */
                margin-right: auto; /* Center tabs in container */
                z-index: 11;

                /* Enable horizontal scrolling */
                overflow-x: scroll; /* Allow horizontal scrolling */
                overflow-y: hidden; /* Prevent vertical scrolling */
                scroll-behavior: smooth; /* Smooth scrolling */
                white-space: nowrap; /* Prevent tabs from wrapping */
                -webkit-overflow-scrolling: touch; /* Smooth scrolling on touch devices */
            }

            /* Ensure the scroll works properly in Firefox */
            div[data-baseweb="tab-list"] {
                scroll-snap-type: x mandatory; /* Ensure proper scrolling behavior */
                scroll-padding: 1rem; /* Add some padding for the scroll */
            }

            div[data-baseweb="tab-border"] {
                height: 0 !important;
            }

            div[data-baseweb="tab-highlight"] {
                background-color: #ff0000;

            }

                        div[data-baseweb="slider"]{
                color: #ff0000;
                }

            div[data-baseweb="slider"] [role="slider"] {
                color: #ff0000;
                fill: #ff0000;
                background-color: #ff0000;
                }

            /* Target the selectbox container */
            div[data-baseweb="select"] {
                background-color: #ffffff !important; /* Light blue background */
                border: 2px solid #b3b3b3 !important; /* Dark blue border */
                border-radius: 10px !important; /* Rounded corners */
                padding: 0.0rem !important; /* Internal padding */
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1) !important; /* Subtle shadow */
            }

                        /* Hover effect for selectboxes */
            div[data-baseweb="select"]:hover {
                border-color: #3d9df3 !important; /* Slightly darker border on hover */
            }

            div[data-baseweb="select"] svg {
                fill: #FF0000;
            }

            div[data-baseweb="select"]:hover svg {
                fill: #00aaf0;
            }

            /* Hover effect for text input */
            div[data-testid="stTextInput"]:hover [data-baseweb="input"] {
                border: 2px solid #3d9df3 !important; /* Slightly darker border on hover */
            }

            div[data-testid="stTooltipIcon"] div[data-testid="stTooltipHoverTarget"] svg {
                fill: none;
            }

            /* Style individual tabs */
            button[data-testid="stTab"] {
                font-size: 0.9rem; /* Slightly smaller text for compactness */
                color: #003366; /* Dark blue text to match your theme */
                background-color: #ffffff; /* Light blue background for tabs */
                border: none; /* Remove all borders */
                border-bottom: 3px solid #99ccff; /* Add a bottom border only */
                border-radius: 0.1rem; /* Slightly rounded edges */
                border-color: #2d9dfa;
                transition: all 0.3s ease-in-out; /* Smooth hover and active effects */
                cursor: pointer;
                padding: 0.3rem 0.3rem !important; /* Increase horizontal padding */

            }

            /* Hover effect for tabs */
            button[data-testid="stTab"]:hover {
                background-color: #92c9fc; /* Slightly darker blue on hover */
                color: #003366; /* Keep text color consistent */
                border-color: #ff0000; /* Highlighted border */
                border-radius: 0.2rem;
                transform: translateY(-1px); /* Subtle lift effect */
            }

            /* Active tab style */
            button[data-testid="stTab"][aria-selected="true"] {
                background-color: #0099ff; /* Bold blue for active tab */
                color: white; /* White text for contrast */
                font-weight: bold;
                border-radius: 0.3rem;
                border-top: 2.5px solid #ff0000; 
                border-bottom: 0 solid #ff0000; 
                transform: translateY(-1px); /* No lift for active tab */
                box-shadow: 0 2px 6px rgba(0, 0, 0, 0.2); /* Shadow for active tab */
            }


            /* Optional: Adjust the header for better alignment */
            header[data-testid="stHeader"] {
                padding-bottom: 10 !important; /* Remove any extra spacing at the bottom */
                padding-top: 10 !important;
                z-index: 10;
                background: transparent;
                height: 0;
            }

            /* Target the st.image container */
            div[data-testid="stImage"] img {
                transition: transform 0.3s ease-in-out; /* Smooth scaling effect */
                cursor: pointer; /* Change cursor to pointer */
            }

            /* On hover, scale the image */
            div[data-testid="stImage"] img:hover {
                transform: scale(1.1); /* Make the image 10% larger */
                z-index: 10;
                box-shadow: 0 4px 8px rgba(0, 0, 0, 0.4); /* Add a shadow effect */
                border: none; /* Remove all borders */
                background-color: #ffffff; /* Remove all borders */
                border-radius: 5px; /* Optional rounded corners */
            }

            div[data-testid="stPlotlyChart"] {
                transition: transform 0.3s ease; /* Smooth scaling transition */
            }
            div[data-testid="stPlotlyChart"]:hover {
                transform: scale(1.025); /* Scale up the chart on hover */
                z-index: 10; /* Ensure it appears above other elements */
                box-shadow: 0 4px 8px rgba(0, 0, 0, 0.4); /* Add a shadow effect */
                border: none; /* Remove all borders */
                background-color: #ffffff; /* Remove all borders */
                border-radius: 5px; /* Optional rounded corners */
            }

                /* Target the expander container */
            div[data-testid="stExpander"] {
                background-color: #ffffff !important; /* Custom light pink background color */
                border-radius: 8px !important; /* Add rounded corners */
                box-shadow: 0 4px 6px rgba(0, 0, 0, 0.2) !important; /* Subtle shadow for depth */
                padding: 0 !important; /* Add internal padding */   
            }

            div[data-testid="stExpander"] span {
                font-size: 16px; /* Make text smaller */
                justify-content: center;

            }

             div[data-testid="stExpander"] summary  [data-testid="stMarkdownContainer"] span {
                font-size: 18px; /* Make text smaller */
                justify-content: center;
                font-weight: bold;
            }

            div[data-testid="stExpander"] svg {
                fill: #ff0000;
            }

            div[data-testid="stSliderThumbValue"] {
                color: #ff0000;
                fill: #ff0000;
            }



            /* Target the st.dataframe container */
            div[data-testid="stDataFrame"] {
                box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2) !important; /* Subtle shadow */
            }


            div[data-testid="stForm"]  {
                box-shadow: 0 3px 6px rgba(0, 0, 0, 0.2) !important; /* Subtle shadow */
            }

            div[data-testid="stCheckbox"]:hover  {
                opacity: 0.5;
            }

            div[data-testid="stWidgetLabel"]:hover svg {
                fill: #FF0000;
            }

            div[data-testid="stMultiSelect"] span[data-baseweb="tag"] svg{
                fill: #ffffff;
            }

            div[data-testid="stMultiSelect"]:hover span[data-baseweb="tag"] svg{
                fill: #ff0000;
            }

            div[data-testid="stMultiSelect"] svg[data-baseweb="icon"] {
                fill: #00aaf0;
            }

            div[data-testid="stMultiSelect"]:hover svg[data-baseweb="icon"] {
                fill: #ff0000;
            }

            div[data-testid="stHeadingWithActionElements"]:hover svg {
                display: none;
            }

            /* Target the Color Picker container */
            div[data-testid="stColorPicker"] {
                transition: transform 0.3s ease-in-out; /* Smooth scaling effect */
            }

            div[data-testid="stColorPicker"]  p{
                font-weight: bold !important;
                color: #0068c9;
            }

            /* Scale up the Color Picker on hover */
            div[data-testid="stColorPickerBlock"]:hover {
                transform: scale(1.25); /* Increase size by 10% */
                z-index: 10; /* Ensure it appears above other elements */
                border-radius: 8px; /* Optional: Add rounded corners */
            }

            /* Target the num inp. container */
            div[data-testid="stNumberInputContainer"] {
                background-color: #ffffff !important; /* Light blue background */
                border: 2px solid #b3b3b3 !important; /* Dark blue border */
                border-radius: 10px !important; /* Rounded corners */
                padding: 0.0rem !important; /* Internal padding */
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1) !important; /* Subtle shadow */
            }


            /* Target the text container */
            div[data-testid="stTextInputRootElement"] {
                background-color: #ffffff !important; /* Light blue background */
                border: 2px solid #b3b3b3 !important; /* Dark blue border */
                border-radius: 10px !important; /* Rounded corners */
                padding: 0.0rem !important; /* Internal padding */
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1) !important; /* Subtle shadow */
            }



            /* Hover effect for number input */
            div[data-testid="stNumberInput"]:hover [data-testid="stNumberInputContainer"] {
                border: 2px solid #3d9df3 !important; /* Slightly darker border on hover */
            }


            /* Close sidebar*/
            div[data-testid="stSidebarCollapsedControl"]{
                margin: -1rem -.75rem;
            }

            /* Close sidebar button color*/
            div[data-testid="stSidebarCollapsedControl"] svg {
                fill: #ff0000;
            }

            /* Hamberger menu*/
            div[data-testid="stToolbar"]{
                margin: -0.5rem -0.5rem;
                color: rgb(255, 0, 0);
            }
            
            div[data-testid="stSpinner"]{
                z-index: 10;
            }

            div[data-testid="stSidebarHeader"] img[data-testid="stLogo"] {
                height: 100px !important; /* Fixed height */
                width: 100% !important; /* Adjust width automatically */
                max-width: 350px !important; /* Restrict the maximum width */
                object-fit: contain !important; /* Proper scaling */
                margin: 0 auto !important; /* Center */
                display: block !important;
                padding: 0 !important;
                visibility: visible !important;
                opacity: 1 !important;
            }
            
            /* Close sidebar*/
            div[data-testid="stSidebarCollapsedControl"] [data-testid="stLogo"]{
                height: 0 !important;
                width: 0 !important;
            }
            
            /* Close sidebar*/
            div[data-testid="stSidebarCollapsedControl"] [data-testid="stLogoLink"]{
                display: none !important;
            }


            div[data-testid="stSidebarCollapseButton"] {
               position: absolute !important; /* Make the position relative to its nearest positioned ancestor */
                top: 1px !important; /* Adjust the distance from the top */
                right: 1px !important; /* Adjust the distance from the right */
                z-index: 1000 !important; /* Ensure it stays on top */
                display: flex !important; /* Ensure it displays as a flexbox */
                align-items: center !important; /* Center the content vertically */
                justify-content: center !important; /* Center the content horizontally */
                color: #ff0000 !important; /* Customize button color */
                width: 1.25rem;
                height: 1.25rem;
                background-color: transparent !important; /* Optional background color */
                border-radius: 5px !important; /* Add rounded corners */
                box-shadow: 0px 1px 2px rgba(0, 0, 0, 0.1) !important; /* Optional shadow effect */
            }

            /* Hover effect for the button */
            div[data-testid="stSidebarCollapseButton"]:hover {
                background-color: #5a93d1 !important; /* Change background on hover */
                color: #ffffff !important; /* Change icon color on hover */
                transform: scale(1.1); /* Slightly enlarge on hover */
                transition: all 0.3s ease-in-out; /* Smooth transition */
                padding: 0 !important; /* Add padding */
                border-radius: 5px !important; /* Add rounded corners */
                margin: 0 !important; /* Margin of the border */
            }

            /* Header divider */
            hr[data-testid="stHeadingDivider"] {
                margin-top: -0.5rem !important;
                margin-bottom: -0.5rem !important;
                border: none;
                border-radius: 1px;
                background-color: #3d9df3;
                height: 0.12rem !important;
            }


        <style>
    """
    st.markdown(custom_style, unsafe_allow_html=True)

    unique_elements_css = """
    <style>
    .st-key-main_img_NE div{
        display: flex;
        justify-content: center;
        align-items: center;
        height: 100%;
    }

    .st-key-subtype_and_or div{
        justify-content: center;
    }

    .st-key-subtype_and_or div svg{
        display:none;
    }

    .st-key-hm_drawer div{
        display: flex;
        justify-content: center;
        align-items: center;
        height: 100%;
    }

    .st-key-hm_drawer button{
        width: 110% !important;
    }

    .st-key-subtype_and_or_tcga div{
        justify-content: center;
    }

    .st-key-subtype_and_or_tcga_sggc div{
        justify-content: center;
    }

    .st-key-subtype_and_or_tcga_sgsg div{
        justify-content: center;
    }

    .st-key-subtype_and_or_sgsg div{
        justify-content: center;
    }

    .st-key-subtype_and_or_sggene_comp div{
        justify-content: center;
    }

    .st-key-subtype_and_or_gf div{
        justify-content: center;
    }

    .st-key-subtype_and_or_gf_tcga div{
        justify-content: center;
    }


    .st-key-main_img_NE img{
        margin: 0;
        max-width: 100%;
        height: auto;
    }

    .st-key-main_img_NE button[data-testid="stBaseButton-elementToolbar"]{
        display: none;
    }

    /* Unique elements CSS */
    .st-key-z_score_toggle_heatmap_heatmap{
        margin-top: 0.65rem !important; /* Add space above to move it lower */
        position: relative; /* Ensure relative positioning for fine-tuning */
        top: 0; /* Move it 20px lower within the parent container */
    }

    .st-key-z_score_toggle_clustering{
        margin-top: 0.65rem !important; /* Add space above to move it lower */
        position: relative; /* Ensure relative positioning for fine-tuning */
        top: 0; /* Move it 20px lower within the parent container */
    }

    <style>
    """

    st.markdown(unique_elements_css, unsafe_allow_html=True)

    custom_script = """
        <style>
            /* Center the running widget */
            div[class*="StatusWidget"] {
                background: rgba(0, 165, 255, 0.5);
                border-radius: 15px;
                height: 175px;
                max-width: 500px;
                width: 350px;
                position: fixed;
                top: 50%;
                left: 50%;
                transform: translate(-50%, -50%);
                display: flex;
                justify-content: center;
                align-items: center;
                opacity: 0.85;
                z-index: 99999;

            }

            /* Change the font style of the text inside the widget */
            div[class*="StatusWidget"] span, 
            div[class*="StatusWidget"] label {
                font-family: Arial, sans-serif !important; /* Use Arial font */
                font-weight: bold !important; /* Make the font bold */
                color: black !important; /* Set the text color to black */
                font-size: 24px !important; /* Optional: Adjust font size */
            }

            div[class*="StatusWidget"] span {
                color: #000000;
            }

            div[class*="StatusWidget"] button {
                display: none;
            }
        </style>
    """
    st.markdown(custom_script, unsafe_allow_html=True)

    # Application logo
    # st.sidebar.image("style_items/surv_sig_logo.svg", use_container_width=True, clamp=True)
    st.logo("style_items/surv_sig_logo.svg", size="large",
            link="https://www.hcemm.eu/teams/genomic-instability-and-cancer/cancer-genomics-and-epigenetics-core-group/",
            icon_image="style_items/favicon.svg")




# Warning message
def warning_tab():
    return st.info(
        ":blue[⚠️ **Application Performance Optimization**⚠️]\n\n"
        ":blue[This application is divided into multiple tabs. If you are not actively using a tab, "
        "you can :red[🚫 turn it off 🚫] using the :blue[toggle] at the top of each tab.]\n\n"
        ":blue[📃 Keeping multiple tabs running simultaneously will ] :red[️📈 increase the overall 🏃‍ running time! ⏱️] \n\n"
        ":blue[✔️ Please use only the tabs you need and turn others off to ensure optimal performance. ⚡ ]"
    )

def display_free_use_message():
    """
    Display an info message styled similarly to the provided example.
    """
    return st.info(
        ":blue[💡 **Free-to-Use Web Application** 💡]\n\n"
        ":blue[This application is completely free to use and no registration is required. Feel free to explore its features "
        "and utilize it for your personal or research purposes.]\n\n"
        ":blue[🔒 Your data remains secure, and the application is designed to provide an optimized experience.]\n\n"
        ":blue[✔️ Thank you for using our application! 😊]"
    )


def landing_animation():
    landing_animation_css = """
    <!-- Overlay container that covers the whole screen -->
    <div class="intro-overlay">
      <div class="logo-container">
        <!-- Your logo -->
    <?xml version="1.0" encoding="UTF-8"?><svg id="Layer_1" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 777.13 278.76"><defs><style>.cls-1{fill:#ff3131;}.cls-2{fill:#009ee2;}.cls-3{fill:none;stroke:#009ee2;stroke-miterlimit:10;stroke-width:2px;}</style></defs><path class="cls-1" d="M44.27,134.4h-20.96c-2.06,0-3.98-1.11-5-2.89l-10.48-18.16c-1.03-1.78-1.03-4,0-5.78l10.48-18.16c1.03-1.78,2.95-2.89,5-2.89h20.96c2.06,0,3.98,1.11,5.01,2.89l10.48,18.16c1.03,1.78,1.03,4,0,5.78l-10.48,18.16c-1.03,1.78-2.95,2.89-5.01,2.89"/><path class="cls-2" d="M90.62,83.03h-34.96c-3.43,0-6.63-1.85-8.34-4.82l-17.48-30.27c-1.72-2.97-1.72-6.67,0-9.64l17.48-30.27c1.72-2.97,4.91-4.82,8.34-4.82h34.96c3.43,0,6.63,1.85,8.35,4.82l17.48,30.27c1.72,2.97,1.72,6.67,0,9.64l-17.48,30.27c-1.72,2.97-4.91,4.82-8.35,4.82"/><path class="cls-1" d="M147.2,102.88c7.74,3.96,19.63,7.93,31.9,7.93,13.21,0,20.2-5.48,20.2-13.78s-6.04-12.46-21.33-17.93c-21.14-7.37-34.92-19.07-34.92-37.57,0-21.71,18.12-38.32,48.14-38.32,14.35,0,24.92,3.02,32.47,6.42l-6.42,23.22c-5.1-2.45-14.16-6.04-26.62-6.04s-18.5,5.66-18.5,12.27c0,8.12,7.17,11.7,23.6,17.93,22.46,8.31,33.03,20.01,33.03,37.94,0,21.33-16.42,39.45-51.34,39.45-14.53,0-28.88-3.77-36.05-7.74l5.85-23.79Z"/><path class="cls-1" d="M330.12,102.88c0,12.08.38,21.9.75,29.64h-24.92l-1.32-13.03h-.56c-3.59,5.66-12.27,15.1-28.88,15.1-18.69,0-32.47-11.7-32.47-40.21v-54.18h28.88v49.64c0,13.41,4.34,21.52,14.35,21.52,7.93,0,12.46-5.47,14.35-10,.75-1.7.94-3.96.94-6.23v-54.93h28.89v62.67Z"/><path class="cls-1" d="M349.04,70.6c0-13.59-.38-22.46-.76-30.39h24.73l.94,16.99h.76c4.72-13.4,16.04-19.06,24.91-19.06,2.65,0,3.97,0,6.04.38v26.99c-2.07-.38-4.53-.75-7.74-.75-10.57,0-17.75,5.66-19.63,14.53-.38,1.89-.57,4.15-.57,6.42v46.81h-28.69v-61.92Z"/><path class="cls-1" d="M449.58,40.21l12.46,42.85c2.26,7.74,3.96,15.1,5.28,22.46h.56c1.51-7.55,3.02-14.53,5.1-22.46l11.89-42.85h30.2l-34.35,92.31h-28.69l-33.6-92.31h31.15Z"/><path class="cls-2" d="M530.17,102.88c7.74,3.96,19.63,7.93,31.9,7.93,13.21,0,20.2-5.48,20.2-13.78s-6.04-12.46-21.33-17.93c-21.14-7.37-34.92-19.07-34.92-37.57,0-21.71,18.12-38.32,48.14-38.32,14.34,0,24.91,3.02,32.46,6.42l-6.42,23.22c-5.1-2.45-14.16-6.04-26.62-6.04s-18.5,5.66-18.5,12.27c0,8.12,7.17,11.7,23.59,17.93,22.46,8.31,33.04,20.01,33.04,37.94,0,21.33-16.42,39.45-51.35,39.45-14.53,0-28.88-3.77-36.05-7.74l5.85-23.79Z"/><rect class="cls-2" x="630.21" y="40.21" width="28.69" height="92.31"/><rect class="cls-2" x="630.21" width="28.69" height="27"/><path class="cls-2" d="M770.08,40.21c-.38,5.85-.76,13.59-.76,27.37v51.53c0,17.74-3.59,32.28-13.97,41.53-10.19,8.68-23.97,11.32-37.56,11.32-12.08,0-24.91-2.45-33.22-7.17l5.66-21.71c5.85,3.4,16.05,6.99,26.81,6.99,13.59,0,23.97-7.36,23.97-24.35v-6.04h-.38c-5.47,7.74-14.35,12.08-24.92,12.08-22.84,0-39.07-18.5-39.07-45.11,0-29.64,19.25-48.51,41.91-48.51,12.65,0,20.58,5.47,25.29,13.02h.38l.94-10.95h24.92ZM740.63,77.77c0-1.89-.19-3.77-.56-5.28-2.08-7.55-7.55-12.65-15.48-12.65-10.38,0-18.88,9.44-18.88,26.24,0,13.78,6.8,24.54,18.88,24.54,7.36,0,13.21-4.91,15.1-11.7.76-2.08.94-5.1.94-7.55v-13.59Z"/><path class="cls-2" d="M286.73,204.77l-30.74-17.74c-1.72-.99-3.64-1.49-5.56-1.49s-3.84.5-5.56,1.49l-30.74,17.74c-3.44,1.99-5.56,5.66-5.56,9.63v35.5c0,3.97,2.12,7.64,5.56,9.63l30.74,17.75c1.72.99,3.64,1.49,5.56,1.49s3.84-.5,5.56-1.49l30.74-17.75c3.44-1.99,5.56-5.66,5.56-9.63v-35.5c0-3.97-2.12-7.64-5.56-9.63M289.6,214.4v3.77c-.43-.67-.91-1.3-1.48-1.88l-18.02-18.02,15.29,8.83c2.6,1.5,4.21,4.3,4.21,7.3M274.44,262.12l-4.59,1.23c1.22-.9,2.27-2.04,3.06-3.4l12.09-20.94-4.6,17.14c-.77,2.9-3.06,5.18-5.96,5.96M236.48,267.93l-3.36-3.36c1.39.61,2.9.95,4.47.95h24.18l-17.15,4.6c-.71.19-1.45.29-2.18.29-2.25,0-4.37-.88-5.96-2.47M226.42,202.17l4.6-1.23c-1.22.9-2.27,2.04-3.06,3.4l-12.09,20.94,4.59-17.14c.78-2.9,3.06-5.18,5.96-5.96M264.38,196.36l3.36,3.36c-1.39-.61-2.9-.95-4.47-.95h-24.18l17.15-4.59c.71-.19,1.44-.29,2.18-.29,2.25,0,4.37.88,5.96,2.47M263.28,261.47h-25.69c-2.52,0-4.87-1.36-6.13-3.54l-12.84-22.25c-1.26-2.19-1.26-4.9,0-7.08l12.84-22.25c1.26-2.18,3.61-3.54,6.13-3.54h25.69c2.52,0,4.87,1.36,6.13,3.54l12.85,22.25c1.26,2.19,1.26,4.9,0,7.08l-12.85,22.25c-1.26,2.18-3.61,3.54-6.13,3.54M213.7,233.36c.17,1.51.63,2.99,1.41,4.35l12.09,20.94-12.55-12.55c-2.12-2.12-2.96-5.24-2.18-8.14l1.23-4.6ZM287.16,230.93c-.17-1.51-.63-2.99-1.41-4.35l-12.09-20.94,12.55,12.55c2.13,2.12,2.96,5.24,2.18,8.14l-1.23,4.6ZM246.22,189.35c1.28-.74,2.74-1.13,4.22-1.13s2.93.39,4.21,1.13l3.27,1.89c-.79.04-1.58.14-2.37.35l-24.61,6.59,15.29-8.83ZM215.48,207.1l3.27-1.88c-.37.7-.68,1.44-.89,2.22l-6.6,24.61v-17.65c0-3,1.62-5.8,4.22-7.3M211.27,249.89v-3.77c.43.67.91,1.3,1.49,1.88l18.02,18.02-15.29-8.83c-2.6-1.5-4.22-4.3-4.22-7.3M254.65,274.94c-1.28.74-2.74,1.13-4.21,1.13s-2.94-.39-4.22-1.13l-3.26-1.89c.79-.03,1.58-.14,2.37-.35l24.61-6.59-15.29,8.83ZM285.38,257.19l-3.27,1.89c.36-.7.67-1.44.88-2.22l6.6-24.61v17.65c0,3-1.61,5.8-4.21,7.3"/><polygon class="cls-2" points="313.22 196.17 323.21 196.17 323.21 217.16 343.21 217.16 343.21 196.17 353.21 196.17 353.21 246.01 343.21 246.01 343.21 225.02 323.21 225.02 323.21 246.01 313.22 246.01 313.22 196.17"/><polygon class="cls-2" points="409.03 196.17 438.45 196.17 438.45 204.02 419.03 204.02 419.03 217.16 436.33 217.16 436.33 225.02 419.03 225.02 419.03 238.15 438.59 238.15 438.59 246.01 409.03 246.01 409.03 196.17"/><polygon class="cls-2" points="447.8 196.31 463.8 196.31 475.66 234.09 475.79 234.09 487.65 196.31 503.65 196.31 503.65 246.15 494.07 246.15 494.07 205.02 493.94 205.02 480.51 246.15 470.94 246.15 457.51 205.02 457.37 205.17 457.37 246.15 447.8 246.15 447.8 196.31"/><polygon class="cls-2" points="512.71 196.17 528.7 196.17 540.56 233.94 540.7 233.94 552.55 196.17 568.55 196.17 568.55 246.01 558.98 246.01 558.98 204.88 558.84 204.88 545.41 246.01 535.84 246.01 522.42 204.88 522.28 205.02 522.28 246.01 512.71 246.01 512.71 196.17"/><path class="cls-2" d="M399.82,237.1c-2.44,1.13-6.82,1.91-10.57,1.91-10.64,0-16.93-7.64-16.93-17.78s6.43-18.07,16.71-18.07c3.57,0,7,.5,10.78,2.64v-.04s0-8.64,0-8.64c-3.49-1.24-7.49-1.82-10.57-1.82-16.93,0-27.35,9.57-27.35,26.63s11.21,24.92,27.35,24.92c3.52,0,7.15-.63,10.57-1.4v-8.37Z"/><polygon class="cls-2" points="327.52 259.46 321.49 259.46 321.49 274.33 319.25 274.33 319.25 259.46 313.22 259.46 313.22 257.46 327.52 257.46 327.52 259.46"/><path class="cls-2" d="M335.25,264h-.11c-.32-.08-.62-.13-.92-.16s-.65-.05-1.06-.05c-.66,0-1.29.14-1.9.44-.61.29-1.2.67-1.77,1.13v8.99h-2.13v-12.66h2.13v1.87c.85-.68,1.59-1.16,2.24-1.44.65-.28,1.31-.42,1.98-.42.37,0,.64,0,.81.03.16.02.41.06.75.11v2.19Z"/><path class="cls-2" d="M346.68,274.34h-2.12v-1.35c-.19.13-.44.31-.77.54-.32.23-.63.41-.93.55-.36.17-.76.32-1.23.44-.46.12-1,.17-1.62.17-1.14,0-2.11-.38-2.9-1.13-.79-.75-1.19-1.72-1.19-2.89,0-.96.2-1.74.62-2.33.41-.59,1-1.06,1.76-1.4.77-.34,1.7-.57,2.78-.69,1.08-.12,2.24-.21,3.48-.27v-.33c0-.49-.09-.88-.25-1.2-.17-.32-.41-.57-.73-.75-.3-.17-.67-.29-1.09-.35-.42-.06-.86-.09-1.32-.09-.56,0-1.18.07-1.87.22-.69.15-1.4.36-2.13.64h-.11v-2.16c.42-.11,1.02-.24,1.8-.37.79-.14,1.56-.2,2.32-.2.89,0,1.67.07,2.33.22.66.15,1.23.4,1.72.75.48.35.84.8,1.09,1.35.25.55.38,1.24.38,2.05v8.59ZM344.56,271.22v-3.52c-.65.04-1.42.09-2.3.17-.88.08-1.57.19-2.09.33-.61.17-1.11.44-1.49.81-.38.37-.56.87-.56,1.51,0,.73.22,1.27.66,1.64.44.37,1.11.55,2.01.55.75,0,1.43-.14,2.05-.44.62-.29,1.19-.64,1.72-1.05"/><path class="cls-2" d="M361.36,274.33h-2.13v-7.21c0-.58-.04-1.13-.1-1.64-.07-.51-.19-.91-.37-1.2-.19-.32-.46-.55-.82-.71-.36-.16-.81-.23-1.38-.23s-1.19.14-1.82.43c-.64.28-1.24.65-1.82,1.1v9.45h-2.13v-12.66h2.13v1.4c.66-.55,1.35-.98,2.06-1.29.71-.31,1.44-.46,2.19-.46,1.37,0,2.41.41,3.13,1.23.72.82,1.08,2.01,1.08,3.56v8.21Z"/><path class="cls-2" d="M374.51,270.69c0,1.16-.48,2.1-1.43,2.85-.95.74-2.26,1.11-3.91,1.11-.94,0-1.8-.11-2.58-.33-.78-.22-1.44-.47-1.97-.73v-2.39h.11c.67.5,1.42.91,2.24,1.21.82.3,1.61.45,2.37.45.94,0,1.67-.15,2.2-.45.53-.3.8-.78.8-1.43,0-.5-.15-.88-.43-1.13-.28-.26-.84-.47-1.65-.66-.3-.07-.7-.15-1.18-.24s-.93-.19-1.33-.3c-1.11-.3-1.9-.73-2.36-1.3-.47-.57-.7-1.27-.7-2.1,0-.52.11-1.01.32-1.47.22-.46.54-.87.98-1.24.42-.36.96-.64,1.62-.84.65-.21,1.38-.31,2.19-.31.75,0,1.52.09,2.29.28.78.18,1.42.41,1.93.67v2.28h-.11c-.54-.4-1.21-.74-1.98-1.01-.78-.28-1.54-.41-2.29-.41s-1.43.15-1.97.45c-.53.3-.8.74-.8,1.33,0,.52.16.91.49,1.18.32.26.83.48,1.54.65.39.09.83.18,1.32.27.49.09.89.17,1.22.25.99.23,1.75.62,2.29,1.17.54.56.8,1.3.8,2.22"/><rect class="cls-2" x="377.58" y="256.7" width="2.13" height="17.63"/><path class="cls-2" d="M393.76,274.34h-2.12v-1.35c-.19.13-.44.31-.76.54-.32.23-.63.41-.93.55-.36.17-.76.32-1.23.44-.46.12-1,.17-1.62.17-1.14,0-2.11-.38-2.9-1.13-.79-.75-1.19-1.72-1.19-2.89,0-.96.2-1.74.62-2.33.41-.59,1-1.06,1.76-1.4.77-.34,1.7-.57,2.78-.69,1.08-.12,2.24-.21,3.48-.27v-.33c0-.49-.09-.88-.25-1.2-.17-.32-.41-.57-.73-.75-.3-.17-.67-.29-1.09-.35-.42-.06-.86-.09-1.32-.09-.56,0-1.18.07-1.87.22-.69.15-1.4.36-2.13.64h-.11v-2.16c.42-.11,1.02-.24,1.8-.37.79-.14,1.56-.2,2.32-.2.89,0,1.67.07,2.33.22.66.15,1.23.4,1.72.75.48.35.84.8,1.09,1.35.25.55.37,1.24.37,2.05v8.59ZM391.64,271.22v-3.52c-.65.04-1.42.09-2.29.17-.88.08-1.58.19-2.09.33-.61.17-1.11.44-1.49.81-.38.37-.56.87-.56,1.51,0,.73.22,1.27.66,1.64.44.37,1.11.55,2.01.55.75,0,1.43-.14,2.05-.44.62-.29,1.19-.64,1.72-1.05"/><path class="cls-2" d="M404.43,274.22c-.4.1-.84.19-1.31.26-.47.07-.89.1-1.26.1-1.29,0-2.27-.35-2.95-1.04-.67-.69-1.01-1.81-1.01-3.34v-6.73h-1.44v-1.79h1.44v-3.64h2.13v3.64h4.4v1.79h-4.4v5.77c0,.67.01,1.19.05,1.56.03.37.14.72.32,1.05.17.3.39.52.69.66.29.14.73.21,1.33.21.35,0,.71-.05,1.09-.15.38-.1.65-.19.82-.25h.11v1.91Z"/><path class="cls-2" d="M409.29,259.56h-2.4v-2.21h2.4v2.21ZM409.16,274.33h-2.13v-12.66h2.13v12.66Z"/><path class="cls-2" d="M424.16,268.01c0,2.06-.53,3.69-1.59,4.88-1.06,1.19-2.47,1.79-4.25,1.79s-3.21-.6-4.27-1.79c-1.05-1.19-1.58-2.82-1.58-4.88s.53-3.69,1.58-4.89c1.05-1.2,2.48-1.8,4.27-1.8s3.19.6,4.25,1.8c1.06,1.2,1.59,2.83,1.59,4.89M421.96,268.01c0-1.64-.32-2.86-.96-3.65-.64-.8-1.53-1.2-2.67-1.2s-2.05.4-2.69,1.2c-.64.8-.96,2.01-.96,3.65s.32,2.79.96,3.61c.64.82,1.54,1.23,2.69,1.23s2.02-.41,2.67-1.22c.64-.81.97-2.02.97-3.62"/><path class="cls-2" d="M438.04,274.33h-2.13v-7.21c0-.58-.04-1.13-.1-1.64-.07-.51-.19-.91-.38-1.2-.19-.32-.46-.55-.81-.71-.36-.16-.81-.23-1.38-.23s-1.19.14-1.82.43c-.64.28-1.24.65-1.82,1.1v9.45h-2.13v-12.66h2.13v1.4c.66-.55,1.35-.98,2.06-1.29.71-.31,1.44-.46,2.19-.46,1.37,0,2.41.41,3.13,1.23.72.82,1.08,2.01,1.08,3.56v8.21Z"/><path class="cls-2" d="M451.98,274.34h-2.12v-1.35c-.19.13-.44.31-.77.54-.32.23-.63.41-.93.55-.36.17-.76.32-1.22.44-.46.12-1,.17-1.62.17-1.14,0-2.11-.38-2.9-1.13-.79-.75-1.19-1.72-1.19-2.89,0-.96.2-1.74.61-2.33.41-.59,1-1.06,1.76-1.4.77-.34,1.69-.57,2.78-.69,1.08-.12,2.24-.21,3.48-.27v-.33c0-.49-.09-.88-.25-1.2-.17-.32-.41-.57-.73-.75-.3-.17-.67-.29-1.09-.35-.42-.06-.86-.09-1.32-.09-.56,0-1.18.07-1.87.22-.69.15-1.4.36-2.13.64h-.11v-2.16c.42-.11,1.02-.24,1.8-.37.79-.14,1.56-.2,2.32-.2.89,0,1.67.07,2.33.22.66.15,1.23.4,1.72.75.47.35.84.8,1.09,1.35.25.55.38,1.24.38,2.05v8.59ZM449.86,271.22v-3.52c-.65.04-1.42.09-2.3.17-.88.08-1.57.19-2.09.33-.61.17-1.11.44-1.49.81-.38.37-.56.87-.56,1.51,0,.73.22,1.27.66,1.64.44.37,1.11.55,2.01.55.75,0,1.43-.14,2.05-.44.62-.29,1.19-.64,1.72-1.05"/><rect class="cls-2" x="456.1" y="256.7" width="2.13" height="17.63"/><polygon class="cls-2" points="485.8 274.33 483.56 274.33 483.56 259.8 478.87 269.69 477.53 269.69 472.87 259.8 472.87 274.33 470.77 274.33 470.77 257.46 473.83 257.46 478.33 266.86 482.68 257.46 485.8 257.46 485.8 274.33"/><path class="cls-2" d="M500.76,268.23h-9.33c0,.78.12,1.46.35,2.03.24.58.56,1.05.96,1.42.39.36.86.63,1.4.82.54.18,1.14.27,1.79.27.86,0,1.73-.17,2.6-.52.87-.34,1.49-.68,1.86-1.01h.11v2.32c-.72.3-1.45.56-2.2.76-.75.2-1.53.31-2.36.31-2.1,0-3.74-.57-4.92-1.71-1.18-1.14-1.77-2.75-1.77-4.84s.56-3.71,1.69-4.93c1.13-1.22,2.62-1.83,4.46-1.83,1.71,0,3.02.5,3.95,1.5.93,1,1.39,2.41,1.39,4.25v1.16ZM498.68,266.6c0-1.12-.29-1.98-.85-2.6-.55-.61-1.4-.92-2.53-.92s-2.05.34-2.73,1.01c-.68.67-1.06,1.51-1.15,2.51h7.25Z"/><path class="cls-2" d="M514.25,274.34h-2.13v-1.33c-.61.53-1.25.94-1.92,1.24-.66.29-1.39.44-2.16.44-1.51,0-2.71-.58-3.6-1.74s-1.33-2.78-1.33-4.84c0-1.07.16-2.03.46-2.87.31-.84.72-1.55,1.24-2.14.51-.57,1.11-1.01,1.8-1.32.68-.3,1.39-.45,2.12-.45.66,0,1.25.07,1.77.21.51.14,1.05.36,1.62.65v-5.48h2.13v17.63ZM512.12,271.22v-7.26c-.57-.26-1.09-.44-1.54-.53-.45-.1-.95-.15-1.49-.15-1.19,0-2.12.42-2.79,1.25-.67.83-1,2.01-1,3.53s.26,2.65.77,3.43c.51.78,1.34,1.17,2.47,1.17.6,0,1.22-.13,1.84-.4.62-.27,1.2-.61,1.73-1.04"/><path class="cls-2" d="M520.73,259.56h-2.4v-2.21h2.4v2.21ZM520.6,274.33h-2.13v-12.66h2.13v12.66Z"/><path class="cls-2" d="M534.17,273.54c-.71.34-1.38.6-2.02.79-.64.19-1.32.28-2.03.28-.91,0-1.75-.13-2.52-.4-.76-.27-1.41-.67-1.96-1.22-.55-.54-.98-1.23-1.28-2.06-.3-.83-.45-1.8-.45-2.91,0-2.07.57-3.69,1.71-4.87,1.13-1.18,2.64-1.77,4.5-1.77.72,0,1.44.1,2.14.31.7.2,1.34.45,1.92.75v2.37h-.11c-.65-.51-1.32-.9-2.01-1.17-.69-.27-1.37-.41-2.02-.41-1.21,0-2.16.41-2.86,1.22-.7.81-1.05,2-1.05,3.57s.34,2.7,1.03,3.52c.68.82,1.64,1.23,2.88,1.23.43,0,.87-.06,1.31-.17.45-.11.85-.26,1.2-.44.31-.16.6-.33.87-.5.27-.18.49-.33.64-.46h.11v2.35Z"/><path class="cls-2" d="M539.19,259.56h-2.4v-2.21h2.4v2.21ZM539.06,274.33h-2.13v-12.66h2.13v12.66Z"/><path class="cls-2" d="M553.85,274.33h-2.13v-7.21c0-.58-.04-1.13-.1-1.64-.07-.51-.19-.91-.38-1.2-.19-.32-.46-.55-.81-.71-.36-.16-.81-.23-1.38-.23s-1.19.14-1.82.43c-.64.28-1.24.65-1.82,1.1v9.45h-2.13v-12.66h2.13v1.4c.66-.55,1.35-.98,2.06-1.29.71-.31,1.44-.46,2.19-.46,1.37,0,2.41.41,3.13,1.23.72.82,1.08,2.01,1.08,3.56v8.21Z"/><path class="cls-2" d="M568.55,268.23h-9.33c0,.78.12,1.46.35,2.03.24.58.56,1.05.96,1.42.39.36.86.63,1.4.82.54.18,1.14.27,1.79.27.86,0,1.73-.17,2.6-.52.87-.34,1.49-.68,1.86-1.01h.11v2.32c-.72.3-1.45.56-2.2.76-.75.2-1.53.31-2.36.31-2.1,0-3.74-.57-4.92-1.71-1.18-1.14-1.77-2.75-1.77-4.84s.56-3.71,1.69-4.93c1.13-1.22,2.62-1.83,4.46-1.83,1.71,0,3.02.5,3.95,1.5.93,1,1.39,2.41,1.39,4.25v1.16ZM566.48,266.6c0-1.12-.29-1.98-.85-2.6-.55-.61-1.4-.92-2.53-.92s-2.05.34-2.73,1.01c-.68.67-1.06,1.51-1.15,2.51h7.25Z"/><line class="cls-3" y1="230.89" x2="159" y2="230.89"/><line class="cls-3" x1="618.12" y1="230.89" x2="777.13" y2="230.89"/></svg>
        <!-- Shine element -->
        <div class="shine"></div>
      </div>
    </div>

    <style>
    /* 1) Full-page overlay */
    .intro-overlay {
      position: fixed;
      top: 25%; 
      left: 0;
      width: 100%;
      height: 100%;
      background-color: #fff; /* Adjust if you want a different color or transparency */
      z-index: 9999;          /* Make sure this is on top of everything */
      display: flex;
      justify-content: center;
      align-items: center;

      /*
        Fades out AFTER the main logo animation + shine are done.
        - spinGrow runs for 3s
        - shine starts at 3s and runs for 1.5s, ending at 4.5s
        This fadeOut animation:
          duration = 1.5s
          delay    = 4.5s
        so total time = 4.5s + 1.5s = 6s before it's fully hidden
      */
      animation: fadeOut 1.5s ease-in 1s forwards;
    }

    /* 2) Container for the logo + shine, so we can position the shine absolutely */
    .logo-container {
      position: relative;
      width: 750px;  /* Adjust as needed */
      height: 750px; /* same as width for a square container */
    }

    /* 3) The logo itself, with spinGrow animation from the center */
    .logo-image {
      width: 300%;
      height: 300%;
      object-fit: contain;
      transform-origin: center center; /* Ensure rotation/grow is from the center */
    }

    /* 4) Shine effect (overlay on top of the logo image) */
    .shine {
      content: "";
      position: absolute;
      top: 0;
      left: -100%;
      width: 100%;
      height: 100%;
      /*
        This gradient creates a diagonal “sweep” highlight
        - Tweak the RGBA values for a more or less intense shine
      */
      background: linear-gradient(
        120deg,
        rgba(255, 255, 255, 0.0) 0%,
        rgba(255, 255, 255, 0.7) 50%,
        rgba(255, 255, 255, 0.0) 100%
      );
      pointer-events: none; /* Don’t block clicks */
      /* Starts after spinGrow finishes (3s), runs for 1.5s */
      animation: shine 1s ease-in-out 1s forwards;
    }


    /* 6) Keyframes for the shine sweeping across */
    @keyframes shine {
      0% {
        left: -100%;
      }
      100% {
        left: 100%;
      }
    }

    /* 7) Keyframes to fade out the entire overlay */
    @keyframes fadeOut {
      0% {
        opacity: 1;
        visibility: visible;
      }
      100% {
        opacity: 0;
        visibility: hidden; /* or display: none; but hidden is usually safer in CSS animations */
      }
    }
    </style>
    """

    st.markdown(landing_animation_css, unsafe_allow_html=True)