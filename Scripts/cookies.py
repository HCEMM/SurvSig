import streamlit as st

def display_cookie_banner():
    """
    Display the cookie banner with the given HTML and CSS.
    """
    st.markdown(
        """
        <style>
        .cookie-bar {
          position: fixed;
          width: 100%;
          bottom: 0;
          top: unset;
          right: 0;
          left: 0;
          height: 50px;
          text-align: center;
          line-height: 30px;
          align-content: center;
          background: #19BCFF;
          color: white;
          font-size: 14px;
          font-family: "Arial", sans-serif;
          font-weight: 100;
          z-index: 99999999;
          transition: transform 0.8s;
          animation: slideIn 0.8s;
          animation-delay: 0.8s;
          margin-bottom: 0;
        }
        @keyframes slideIn {
          0% {
            transform: translateY(-50px);
          }
          100% {
            transform: translateY(0);
          }
        }
        .close-cb {
            display: inline-flex !important;
            -webkit-box-align: center !important;
            align-items: center !important;
            -webkit-box-pack: center !important;
            justify-content: center !important;
            font-weight: 500 !important;
            padding: 0.5rem 1rem !important; 
            border-radius: 0.5rem !important;
            min-height: 2.5rem !important;
            margin: 0px !important;
            color: rgb(0, 153, 255);
            line-height: 1.6 !important;
            font-weight: bold !important;
            width: auto !important;
            user-select: none !important;
            background-color: white !important; 
            border: 2px solid rgb(255 0 0) !important; 
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
            transition: all 0.3s ease-in-out; 
            height: 80% !important;
            cursor: pointer;
            margin-left: 2% !important;
        }
        .close-cb:hover {
            border-color: rgb(163 0 0) !important; /* Softer light blue */
            box-shadow: 0 5px 8px rgba(0, 0, 0, 0.4); /* Slightly stronger shadow */
            transform: translateY(-2px); /* Lift effect */
        }
        .checkbox-cb {
          display: none;
        }
        
        .message + .cookie-bar{
          color: white;
        }
        
        .checkbox-cb:checked + .cookie-bar {
          transform: translateY(100%);
        }
        a {
          color: rgb(92, 242, 255);
        }
        </style>

        <input class="checkbox-cb" id="checkbox-cb" type="checkbox" />
        <div class="cookie-bar">
        <span class="message">We use cookies to enhance your browsing experience. See our <a href="https://www.hcemm.eu/cookies-privacy-policy/"><b>Privacy Policy</b></a>.</span>
        <label for="checkbox-cb" class="close-cb">Accept</label>
        </div>
        """,
        unsafe_allow_html=True,
    )
